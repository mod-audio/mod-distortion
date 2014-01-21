
# run make in each plugins subdirectory
.PHONY:
	all

all:
	cd DS1 && make
	cd guitarix-Overdrive && make

install:
	cd DS1 && make INSTALL_PATH=${LV2_PATH} install
	cd guitarix-Overdrive && make INSTALL_PATH=${LV2_PATH} install

clean:
	cd DS1 && make clean
	cd guitarix-Overdrive && make clean
