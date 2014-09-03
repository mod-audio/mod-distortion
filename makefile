
# run make in each plugins subdirectory

LV2_PATH:=/usr/local/lib/lv2

.PHONY:
	all

all:
	cd ds1 && make
	cd bigmuff && make

install:
	cd ds1 && make INSTALL_PATH=${LV2_PATH} install
	cd bigmuff && make INSTALL_PATH=${LV2_PATH} install

clean:
	cd ds1 && make clean
	cd bigmuff && make clean
