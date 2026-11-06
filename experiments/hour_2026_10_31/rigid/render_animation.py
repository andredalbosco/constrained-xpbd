"""Render the saved modal-transfer experiment; no solver or browser is involved."""
from pathlib import Path
import csv
import hashlib
import json
import math
import subprocess
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "animation"
OUT.mkdir(exist_ok=True)
WIDTH, HEIGHT, FPS = 1920, 1080, 30
MAGNIFICATION = 16
FONT = "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"
BOLD = "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"
fonts = {size: ImageFont.truetype(FONT, size) for size in (18, 20, 22, 24, 28, 38)}
bold = {size: ImageFont.truetype(BOLD, size) for size in (20, 22, 24, 38)}
BACKGROUND, INK, MUTED = "#f6f7f9", "#162a34", "#526675"
MODE_COLORS = ("#e07a39", "#8d5cb6", "#147f91")
METHODS = [
    ("reference", "Refined nonlinear", "reference", "#273b49"),
    ("obj_midpotential_b1", "Midpoint potential", "Old reaction omitted", "#d1495b"),
    ("obj_quarter_b1", "Quarter objective", "Old reaction included", "#00798c"),
    ("canonical_dg", "Canonical discrete", "gradient", "#bc850e"),
]
rows = []
for method, *_ in METHODS:
    with (ROOT / "modal_transfer" / (method + ".csv")).open() as stream:
        rows.append(list(csv.DictReader(stream)))
assert all(len(series) == 1800 for series in rows)
for series in rows:
    for index, row in enumerate(series):
        assert abs(float(row["t"]) - (index + 1) / FPS) < 1e-10


def centered(draw, x, y, text, font, fill=INK):
    box = draw.textbbox((0, 0), text, font=font)
    draw.text((x - (box[2] - box[0]) / 2, y), text, font=font, fill=fill)


def render(index):
    image = Image.new("RGB", (WIDTH, HEIGHT), BACKGROUND)
    draw = ImageDraw.Draw(image)
    draw.text((40, 28), "Constant total energy can hide transfer to the wrong mode", font=bold[38], fill=INK)
    draw.text((40, 87), "Saved triple-pendulum trajectories  |  30 Hz simulation  |  1x playback", font=fonts[24], fill=MUTED)
    draw.text((40, 123), f"Rod angles magnified {MAGNIFICATION}x for visibility. Energy diagnostics use the original simulation.", font=fonts[24], fill=MUTED)
    centered(draw, 1710, 100, f"t = {float(rows[0][index]['t']):05.2f} / 60.00 s", bold[24])
    for panel, ((method, title, subtitle, color), series) in enumerate(zip(METHODS, rows)):
        left = 40 + panel * 470
        right = left + 430
        center = (left + right) / 2
        draw.rounded_rectangle((left, 181, right, 1017), radius=18, fill="white", outline="#dce3e8", width=2)
        centered(draw, center, 200, title, bold[24], color)
        centered(draw, center, 236, subtitle, fonts[22], color)
        row = series[index]
        origin_y, scale = 307, 150
        draw.line((center, origin_y, center, 780), fill="#e1e7eb", width=2)
        old_x = old_y = x = y = 0.0
        points = [(center, origin_y)]
        for link in range(3):
            next_x, next_y = float(row[f"x{2*link}"]), float(row[f"x{2*link+1}"])
            angle = math.atan2(next_x - old_x, -(next_y - old_y)) * MAGNIFICATION
            x += math.sin(angle)
            y -= math.cos(angle)
            old_x, old_y = next_x, next_y
            px, py = center + scale * x, origin_y - scale * y
            assert left + 20 <= px <= right - 20 and 280 <= py <= 780
            points.append((px, py))
        draw.line(points, fill=color, width=7, joint="curve")
        for px, py in points[1:]:
            draw.ellipse((px - 11, py - 11, px + 11, py + 11), fill=color, outline="white", width=2)
        draw.ellipse((center - 8, origin_y - 8, center + 8, origin_y + 8), fill=INK)
        centered(draw, center, 798, "Linear-mode energy / initial high-mode energy", fonts[18])
        for mode in range(3):
            ratio = float(row[f"mode{mode}_energy"])
            yy = 836 + mode * 37
            draw.text((left + 22, yy - 3), f"{mode + 1}", font=bold[20], fill=MODE_COLORS[mode])
            draw.rounded_rectangle((left + 53, yy, left + 295, yy + 18), radius=4, fill="#edf0f3")
            if ratio > 0:
                draw.rounded_rectangle((left + 53, yy, left + 53 + 242 * min(ratio, 1), yy + 18), radius=4, fill=MODE_COLORS[mode])
            draw.text((left + 305, yy - 4), f"{ratio:.2e}", font=fonts[18], fill=INK)
        centered(draw, center, 956, f"Total energy residual: {float(row['relative_energy']):+.2e}", fonts[18])
        centered(draw, center, 984, "Same fixed viewport and display scale", fonts[18], MUTED)
    draw.text((40, 1034), "Initial highest-mode amplitude 0.02; low-mode seed 2e-7. Rod lengths remain one unit in the display.", font=fonts[20], fill=MUTED)
    return image


movie = OUT / "modal_transfer_60s.mp4"
command = ["ffmpeg", "-hide_banner", "-loglevel", "error", "-y", "-f", "rawvideo", "-vcodec", "rawvideo", "-pixel_format", "rgb24", "-video_size", f"{WIDTH}x{HEIGHT}", "-framerate", str(FPS), "-i", "-", "-an", "-c:v", "libx264", "-preset", "fast", "-crf", "20", "-pix_fmt", "yuv420p", "-threads", "2", "-movflags", "+faststart", str(movie)]
selected = {0, 599, 1199, 1799}
stills = []
with subprocess.Popen(command, stdin=subprocess.PIPE) as encoder:
    for index in range(1800):
        frame = render(index)
        if index in selected:
            frame.save(OUT / f"frame_{index+1:04d}.png")
            stills.append(frame)
        encoder.stdin.write(frame.tobytes())
        if (index + 1) % 300 == 0:
            print(f"Rendered {index + 1}/1800 frames", flush=True)
    encoder.stdin.close()
    assert encoder.wait() == 0
montage = Image.new("RGB", (2 * WIDTH, 2 * HEIGHT), BACKGROUND)
for index, frame in enumerate(stills):
    montage.paste(frame, ((index % 2) * WIDTH, (index // 2) * HEIGHT))
montage.save(OUT / "representative_frames.png")
probe = json.loads(subprocess.check_output(["ffprobe", "-v", "error", "-count_frames", "-show_entries", "stream=codec_name,width,height,r_frame_rate,nb_read_frames,pix_fmt,duration", "-of", "json", str(movie)]))
stream = probe["streams"][0]
assert stream["nb_read_frames"] == "1800" and abs(float(stream["duration"]) - 60) < 1e-6
subprocess.run(["ffmpeg", "-hide_banner", "-v", "error", "-i", str(movie), "-f", "null", "-"], check=True)
manifest = {"source": "Recorded modal_transfer CSVs; no new simulation", "angle_magnification": MAGNIFICATION, "playback_speed": 1, "shared_fixed_viewport": True, "representative_simulation_times": [1/30, 20, 40, 60], "ffprobe": probe, "full_decode": "pass", "inputs_sha256": {method: hashlib.sha256((ROOT / "modal_transfer" / (method + ".csv")).read_bytes()).hexdigest() for method, *_ in METHODS}, "video_sha256": hashlib.sha256(movie.read_bytes()).hexdigest()}
(OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
print("Saved", movie, flush=True)
