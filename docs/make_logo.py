import point_cloud_utils as pcu
import numpy as np
from scipy.ndimage import gaussian_filter
import argparse
from PIL import Image
from mayavi import mlab
import matplotlib.pyplot as plt


def pad_image(im, pct_lr, pct_ud):
    pad_lr = int(np.round(im.shape[0] * pct_lr))
    pad_ud = int(np.round(im.shape[1] * pct_ud))
    pad_im = np.full((im.shape[0] + 2*pad_lr, im.shape[1] + 2*pad_ud), 0.0)
    pad_im[pad_lr:pad_lr+im.shape[0], pad_ud:pad_ud+im.shape[1]] = im
    return pad_im


def lloyd_sample_im(im, npts):
    uv = pcu.lloyd_2d(npts)
    idx = np.floor(uv * np.array(im.shape)).astype(np.int32)
    h = []
    for i in range(idx.shape[0]):
        h.append(im[idx[i, 0], idx[i, 1]])
    h = np.array(h)
    return np.concatenate([idx/max(*im.shape), h[:, np.newaxis]], axis=-1)


def main():
    argparser = argparse.ArgumentParser()
    argparser.add_argument("text_img")
    argparser.add_argument("--filter", type=int, default=5.0)
    argparser.add_argument("--max-height", type=float, default=0.01)
    argparser.add_argument("--npts", type=int, default=8192)
    argparser.add_argument("--sf", type=float, default=0.005)
    argparser.add_argument("--cm", type=str, default="magma")
    argparser.add_argument("--flip", action="store_true")
    argparser.add_argument("--trim", action="store_true")
    argparser.add_argument("--floor", type=float, default=0.0)

    args = argparser.parse_args()

    im = args.max_height * (255.0 - np.asarray(Image.open(args.text_img))) / 255.0
    im_pad = pad_image(im, 0.2, 0.05)
    im_filt = gaussian_filter(im_pad, args.filter)
    if args.flip:
        im_filt = im_filt.max() - im_filt

    # plt.imshow(im_filt)
    # plt.colorbar()
    # plt.show()
    pts = lloyd_sample_im(im_filt, args.npts)
    clr = np.copy(pts[:, -1])
    if args.flip:
        pts[:, -1] = pts[:, -1].max() - pts[:, -1]

    if args.trim:
        mask = clr > 1e-7
        pts = pts[mask]
        clr = clr[mask]
    mlab.figure(bgcolor=(1.0, 1.0, 1.0), fgcolor=(0.2, 0.2, 0.2))
    mlab.points3d(pts[:, 0], pts[:, 1], pts[:, 2], clr + args.floor * args.max_height, colormap=args.cm,
                  scale_factor=args.sf, vmin=0.0)
    mlab.show()


if __name__ == "__main__":
    main()


