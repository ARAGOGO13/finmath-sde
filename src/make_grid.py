#!/usr/bin/env python3
"""
Create a PDF grid from input PDFs.

Layouts:
  2x2      – 4 images in 2 columns × 2 rows
  2x2+1    – 5 images: 2+2 rows then 1 centered below
  default  – auto-detect from number of files (≤4 → 2x2, 5 → 2x2+1)

Usage:
  make_grid.py output.pdf [--layout 2x2|2x2+1] file1.pdf [file2.pdf ...]

Requires: Pillow, pdftoppm (poppler)
"""
import sys
import subprocess
import os
import tempfile
from PIL import Image


DPI = 250
GAP = 8  # px between cells


def pdf_to_image(pdf_path: str, tmpdir: str) -> Image.Image:
    prefix = os.path.join(tmpdir, os.path.basename(pdf_path).replace(' ', '_'))
    subprocess.run(
        ['pdftoppm', '-r', str(DPI), '-png', '-singlefile', pdf_path, prefix],
        check=True, capture_output=True
    )
    return Image.open(prefix + '.png').convert('RGB')


def _load_images(files: list, tmpdir: str) -> list:
    imgs = []
    for f in files:
        if os.path.exists(f):
            imgs.append(pdf_to_image(f, tmpdir))
    return imgs


def _uniform_size(imgs: list) -> tuple:
    W = max(img.width  for img in imgs)
    H = max(img.height for img in imgs)
    return W, H


def _resize_all(imgs: list, W: int, H: int) -> list:
    return [img.resize((W, H), Image.LANCZOS) for img in imgs]


def make_grid_2x2(files: list, output: str) -> bool:
    """4 images in a 2×2 grid."""
    with tempfile.TemporaryDirectory() as tmpdir:
        imgs = _load_images(files, tmpdir)
        if not imgs:
            return False
        while len(imgs) < 4:
            imgs.append(Image.new('RGB', imgs[0].size, 'white'))
        W, H = _uniform_size(imgs)
        imgs = _resize_all(imgs[:4], W, H)

        grid = Image.new('RGB', (2*W + GAP, 2*H + GAP), 'white')
        grid.paste(imgs[0], (0,       0))
        grid.paste(imgs[1], (W + GAP, 0))
        grid.paste(imgs[2], (0,       H + GAP))
        grid.paste(imgs[3], (W + GAP, H + GAP))
        grid.save(output, 'PDF', resolution=DPI)
    return True


def make_grid_2x2_plus1(files: list, output: str) -> bool:
    """5 images: two full rows (2+2) then one image centered in a third row."""
    with tempfile.TemporaryDirectory() as tmpdir:
        imgs = _load_images(files, tmpdir)
        if not imgs:
            return False
        while len(imgs) < 5:
            imgs.append(Image.new('RGB', imgs[0].size, 'white'))
        W, H = _uniform_size(imgs)
        imgs = _resize_all(imgs[:5], W, H)

        total_w = 2*W + GAP
        total_h = 3*H + 2*GAP
        grid = Image.new('RGB', (total_w, total_h), 'white')

        # Row 1
        grid.paste(imgs[0], (0,       0))
        grid.paste(imgs[1], (W + GAP, 0))
        # Row 2
        grid.paste(imgs[2], (0,       H + GAP))
        grid.paste(imgs[3], (W + GAP, H + GAP))
        # Row 3 — centered
        x_center = (total_w - W) // 2
        grid.paste(imgs[4], (x_center, 2*(H + GAP)))

        grid.save(output, 'PDF', resolution=DPI)
    return True


def main():
    args = sys.argv[1:]
    if len(args) < 2:
        print(f'Usage: {sys.argv[0]} output.pdf [--layout 2x2|2x2+1] file1.pdf ...', file=sys.stderr)
        sys.exit(1)

    output_file = args[0]
    rest = args[1:]

    layout = None
    if rest and rest[0] == '--layout':
        layout = rest[1]
        input_files = rest[2:]
    else:
        input_files = rest

    # Auto-detect layout
    if layout is None:
        layout = '2x2+1' if len(input_files) == 5 else '2x2'

    if layout == '2x2+1':
        ok = make_grid_2x2_plus1(input_files, output_file)
    else:
        ok = make_grid_2x2(input_files, output_file)

    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()