
# source files directory
SRC_DIR = .

# extension source file
EXT = cpp

# compiler and linker
CC = g++
LD = g++

# flags
CFLAGS = -O3 -funroll-loops -ffast-math -fomit-frame-pointer -fstrength-reduce -fPIC -DPIC $(SSE_CFLAGS) $(SSE3_CFLAGS)
LDFLAGS = -shared -nostartfiles -Wl,-Bsymbolic

ifneq (cat /proc/cpuinfo | grep sse >/dev/null,)
SSE_CFLAGS = -msse -mfpmath=sse
endif

ifneq (cat /proc/cpuinfo | grep ssse3 >/dev/null,)
SSE3_CFLAGS = -msse3
endif

# libs
LIBS = -lc -lm -lrt

# remove command
RM = rm -f

# package name
PACKAGE = $(shell cd .. && basename `pwd` | tr A-Z a-z)

# plugin name
PLUGIN = overdrive
PLUGIN_SO = $(PLUGIN).so

# effect path
EFFECT_PATH = $(PLUGIN).lv2

# installation path
ifndef INSTALL_PATH
INSTALL_PATH = /usr/local/lib/lv2
endif
INSTALLATION_PATH = $(INSTALL_PATH)/$(EFFECT_PATH)


# sources and objects
SRC = $(wildcard $(SRC_DIR)/*.$(EXT))
OBJ = $(SRC:.$(EXT)=.o)

## rules
$(PLUGIN): $(OBJ)
	$(LD) $(LDFLAGS) $(CFLAGS) $(OBJ) -o $(PLUGIN_SO) $(LIBS)

%.o: %.$(EXT)
	$(CC) $(CFLAGS) -DURI_PREFIX='"http://portalmod.com/plugins/guitarix"' -c -o $@ $<

clean:
	$(RM) $(SRC_DIR)/*.o *.so $(PLUGIN)

install:
	mkdir -p $(INSTALLATION_PATH)
	cp $(PLUGIN_SO) $(INSTALLATION_PATH)
	cp $(SRC_DIR)/*.ttl $(INSTALLATION_PATH)
	mkdir -p $(INSTALLATION_PATH)/modgui
	cp modgui/*.html modgui/*.json $(INSTALLATION_PATH)/modgui
