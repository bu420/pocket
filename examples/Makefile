TARGETS = window text cube model minecraft flightsim

SRC = ../src/pok_core.c ../src/pok_math.c ../src/pok_util.c ../src/pok_gfx.c ../src/pok_io.c ../src/pok_platform.c
WARN = -Wall -Wextra -Wno-missing-braces -Wno-unused-parameter -Wno-unused-variable

# Debug/release build options.
# Add MODE=RELEASE to enable optimization.
ifeq ($(MODE), RELEASE)
OPTIONS = -DNDEBUG -O3
else
OPTIONS = -g
endif

$(TARGETS): %: %.c
	gcc $^ $(SRC) $(WARN) $(OPTIONS) -I../inc/ -lgdi32 -o $@

all: $(TARGETS)

.PHONY: clean
clean:
	del *.exe
	del *.bmp
