
# run make in each plugins subdirectory
.PHONY:
	all

all:
	cd ds1 && make
	cd guitarix-Overdrive && make
	cd bigmuff && make
	cd mufffuzz && make

install:
	cd ds1 && make INSTALL_PATH=${LV2_PATH} install
	cd guitarix-Overdrive && make INSTALL_PATH=${LV2_PATH} install
	cd bigmuff && make INSTALL_PATH=${LV2_PATH} install
	cd mufffuzz && make INSTALL_PATH=${LV2_PATH} install

clean:
	cd ds1 && make clean
	cd guitarix-Overdrive && make clean
	cd bigmuff && make clean
	cd mufffuzz && make clean
