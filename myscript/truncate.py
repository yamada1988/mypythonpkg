from PIL import Image, ImageDraw, ImageFont

##
# Ref. site: https://note.nkmk.me/pillow/
#

def add_margin(pil_img, top, right, bottom, left, color):
    width, height = pil_img.size
    new_width = width + right + left
    new_height = height + top + bottom
    result = Image.new(pil_img.mode, (new_width, new_height), color)
    result.paste(pil_img, (left, top))
    return result

def get_concat_tile_resize(im_list_2d, resample=Image.BICUBIC):
    im_list_v = [get_concat_h_multi_resize(im_list_h, resample=resample) for im_list_h in im_list_2d]
    return get_concat_v_multi_resize(im_list_v, resample=resample)

def get_concat_h_multi_resize(im_list, resample=Image.BICUBIC):
    min_height = min(im.height for im in im_list)
    im_list_resize = [im.resize((int(im.width * min_height / im.height), min_height),resample=resample)
                      for im in im_list]
    total_width = sum(im.width for im in im_list_resize)
    dst = Image.new('RGB', (total_width, min_height))
    pos_x = 0
    for im in im_list_resize:
        dst.paste(im, (pos_x, 0))
        pos_x += im.width
    return dst

def get_concat_v_multi_resize(im_list, resample=Image.BICUBIC):
    min_width = min(im.width for im in im_list)
    im_list_resize = [im.resize((min_width, int(im.height * min_width / im.width)),resample=resample)
                      for im in im_list]
    total_height = sum(im.height for im in im_list_resize)
    dst = Image.new('RGB', (min_width, total_height))
    pos_y = 0
    for im in im_list_resize:
        dst.paste(im, (0, pos_y))
        pos_y += im.height
    return dst

dt = 100
frames = 2
timages = []
for j in range(1, frames+1):
    images = []
    for i in range(1, 31):
        im = Image.open('md{0:03d}_{1:d}.eps'.format(i, j))
        im_crop = im.crop((30, 40, 480, 500))
        draw = ImageDraw.Draw(im_crop)
        font = ImageFont.truetype("/Library/Fonts/Arial.ttf", 24)
        draw.text((360, 20), 'sys{0:02d}'.format(i), font = font, fill='#000000')
        images.append(im_crop)

    get_concat_tile_resize([images[0:6],
                            images[6:12],
                            images[12:18],
                            images[18:24],
                            images[24:30]]).save('trancate{0:03d}.eps'.format(j))

    im = Image.open('trancate{0:03d}.eps'.format(j))
    im_new = add_margin(im, 80, 0, 0, 0, (255, 255, 255))
    draw = ImageDraw.Draw(im_new)

    font = ImageFont.truetype("/Library/Fonts/Arial.ttf",54)
    draw.text((50, 10), 't = {0:3d}~{1:3d} ps'.format((j-1)*dt, j*dt), font = font, fill='#000000') 

    font = ImageFont.truetype("/Library/Fonts/Arial.ttf",60)
    draw.text((640, 10), 'MOL NUMBER DensMap'.format((j-1)*dt, j*dt), font = font, fill='#000000')
    im_new.save('trancate{0:03d}.eps'.format(j))

    timages.append(im_new)

timages[0].save('trancate.gif', save_all=True, append_images=timages[1:], optimize=False, duration=200, loop=10)
