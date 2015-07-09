
NAME = mod-bigmuff

# installation path
ifndef INSTALL_PATH
INSTALL_PATH = /usr/local/lib/lv2
endif
COMPLETE_INSTALL_PATH = $(DESTDIR)$(INSTALL_PATH)/$(NAME).lv2

# compiler and linker
CXX ?= g++

# flags
CXXFLAGS += -I. -O3 -ffast-math -Wall -Wextra -c -fPIC -DPIC
LDFLAGS += -shared -lm -fPIC -DPIC

# remove command
RM = rm -f

# sources and objects
SRC = $(wildcard src/*.cpp)
OBJ = $(SRC:.cpp=.o)

# plugin name
PLUGIN = $(NAME).so

$(PLUGIN): $(OBJ)
	$(CXX) $(LDFLAGS) $(OBJ) -o $(PLUGIN)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	$(RM) src/*.o $(PLUGIN)

install:
	mkdir -p $(COMPLETE_INSTALL_PATH)
	cp $(PLUGIN) $(COMPLETE_INSTALL_PATH)
	cp -r ttl/* $(COMPLETE_INSTALL_PATH)
