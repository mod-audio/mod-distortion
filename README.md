mod-distortion
==============

Analog distortion emulation developed by mod team (lv2)

The effects were developed suposing you have a -15dB input signal (measured with digital peak meter)
when you play loud, so its recommended that you ajust your input gain to this level.

We recommend that you use only the stable plugins:
-DS1
-Big Muff Pi

Installation:

	make
	sudo make install
	
You also can install only the plugins you want.
To do this, enter the folder of the effect that you want and type the following comands:

	make
	sudo make install

In this way, the default instalation path is /usr/local/lib/lv2/, and this can be modified passing the variable INSTALL_PATH to make install, e.g.:

	sudo make install INSTALL_PATH=/usr/lib/lv2
