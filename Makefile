
# run make in each plugins subdirectory

.PHONY:
	all

all:
	$(MAKE) -C ds1
	$(MAKE) -C bigmuff

install:
	$(MAKE) -C ds1 install
	$(MAKE) -C bigmuff install

clean:
	$(MAKE) -C ds1 clean
	$(MAKE) -C bigmuff clean
