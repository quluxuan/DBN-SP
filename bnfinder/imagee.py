
from PIL import Image, ImageFilter

#from pylab import *
import numpy as np
import matplotlib.pyplot as plt
#im = Image.open('/home/bubble/Firefox_wallpaper.png')

#w, h = im.size
#print('Original image size: %sx%s' % (w, h))

#im.thumbnail((w//2, h//2))
#print('Resize image to: %sx%s' % (w//2, h//2))
#im2 = im.filter(ImageFilter.BLUR)
#im.save('Firefox_wallpaper.png', 'png')
#im2.show()


x = np.linspace(-np.pi, np.pi, 256, endpoint=True)
C, S = np.cos(x), np.sin(x)
plt.plot(x, C)
plt.plot(x, S)
plt.plot(x, C, color='red', linewidth=2.5, linestyle='-')
plt.plot(x, S, color='yellow', linewidth=2.5, linestyle='-')
plt.show()