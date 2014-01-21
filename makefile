
# run make in each plugins subdirectory
.PHONY:
	all

all:
	cd DS1 && make

install:
	cd DS1 && make INSTALL_PATH=${LV2_PATH} install

clean:
	cd DS1 && make clean
