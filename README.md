# PanoSlice

Upload a dental panoramic X-ray (OPG) and explore an estimated 3D jaw with linked axial, coronal, sagittal and cross-section slices.

Everything runs in the browser. Images are never uploaded, and patient details in DICOM files are not read.

## How it works

- A small U-Net (1.1 M parameters) finds the teeth on the OPG. It was trained on 116 panoramic radiographs and scores a Dice of 0.92 on 20 held-out images.
- The OPG is wrapped around an average adult dental arch. Depth, gums and bone shape come from typical anatomy, not from the image.
- The volume is turned into surface meshes (teeth, bone, gums) in a background thread. Drawing meshes is light work for any graphics chip, so rotation stays smooth on laptops with integrated graphics and on phones. It uses WebGL 2, falls back to WebGL 1, and falls back again to a software renderer when WebGL is switched off.
- The 3D model appears within a few seconds. The U-Net then refines the teeth in the background and the model updates itself.
- Phones and low-memory devices get a lighter model automatically. Add `?full` to the address to force the full model, or `?lite` to force the light one.

## Not a medical device

PanoSlice is for teaching and demos. Don't use it for diagnosis or treatment planning. For real 3D imaging, take a CBCT.

## Credits and licence

- Training images: H. Abdi, S. Kasaei and M. Mehdizadeh, "Automatic segmentation of mandible in panoramic x-ray", J. Med. Imaging 2(4), 2015.
- Tooth masks: S. Helli and A. Hamamcı (2022), via the SerdarHelli/SegmentationOfTeethPanoramicXRayImages dataset on Hugging Face.
- That data is licensed for non-commercial research use, and so are the model weights in `assets/`.

## Deploy

It's a static site with no build step. Keep `index.html` and the `assets/` folder together.

- GitHub Pages: put these files in the repository (or its `docs/` folder), then turn on Pages under Settings → Pages.
- Vercel: import the repo and deploy; `vercel.json` sets the security and cache headers.
- Netlify: `_headers` and `netlify.toml` do the same.
