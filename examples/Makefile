# Targets that only include psr.h
PSR_TARGETS = text
# Targets that include both psr.h and pwa.h
PSR_AND_PWA_TARGETS = window cube model minecraft

WARN = -Wall -Wextra -Wno-missing-braces -Wno-unused-parameter -Wno-unused-variable
OPTIONS = 

# Debug/release build options.
ifeq ($(MODE), RELEASE)
OPTIONS += -DNDEBUG -O3
else
OPTIONS += -g
endif

$(PSR_TARGETS): %: %.c
	gcc $^ ../src/psr.c $(WARN) $(OPTIONS) -I../inc/ -o $@

$(PSR_AND_PWA_TARGETS): %: %.c
	gcc $^ ../src/psr.c ../src/pwa.c $(WARN) $(OPTIONS) -I../inc/ -lgdi32 -o $@

all: $(PSR_TARGETS) $(PSR_AND_PWA_TARGETS)

.PHONY: clean
clean:
	del *.exe
	del *.bmp
