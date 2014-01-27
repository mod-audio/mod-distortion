
NAME = mod-DS1

# installation path
INSTALL_PATH = /usr/local/lib/lv2
COMPLETE_INSTALL_PATH = $(INSTALL_PATH)/$(NAME).lv2

# compiler and linker
CC = g++
LD = g++

# flags
CFLAGS = -I. -O6 -ffast-math -Wall -Wextra -march=native -mtune=native -c -fPIC -DPIC
LDFLAGS = -shared -lc -lm -fPIC -DPIC

# libs
LIBS = -lm

# remove command
RM = rm -f

# sources and objects
SRC = $(wildcard src/*.cpp)
OBJ = $(SRC:.cpp=.o)

# plugin name
PLUGIN = $(NAME).so

$(PLUGIN): $(OBJ)
	$(LD) $(LDFLAGS) $(OBJ) -o $(PLUGIN) $(LIBS)

%.o: %.cpp
	$(CC) $(CFLAGS) -o $@ $<

clean:
	$(RM) src/*.o $(PLUGIN)

install:
	mkdir -p $(COMPLETE_INSTALL_PATH)
	cp $(PLUGIN) $(COMPLETE_INSTALL_PATH)
	cp ttl/*.ttl $(COMPLETE_INSTALL_PATH)
	#mkdir $(COMPLETE_INSTALL_PATH)/modgui
	#cp ttl/modgui/*.html ttl/modgui/*.json $(COMPLETE_INSTALL_PATH)/modgui
