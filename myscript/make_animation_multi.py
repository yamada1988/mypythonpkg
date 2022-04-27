import os
from PIL import Image, ImageDraw, ImageFont

dt = 100
Ntot = 10

font = ImageFont.truetype("/Library/Fonts/Arial.ttf", 24)

for j in range(1, 33):
    images = []
    for i in range(0, Ntot):
        fname = 'md{0:03d}_{1:d}'.format(j, i+1)+'.eps'
        if not os.path.exists(fname):
            break

        im = Image.open(fname)
        draw = ImageDraw.Draw(im)
        draw.text((300, 60), 't = {0:3d}~{1:3d} ps'.format(i*dt, (i+1)*dt), font = font, fill='#000000')
        images.append(im)

    images[0].save('md{0:03d}.gif'.format(j), save_all=True, append_images=images[1:], optimize=False, duration=200, loop=10)
