PSR_TARGETS = text
PSR_AND_PWA_TARGETS = window model minecraft

WARN = -Wall -Wextra -Wno-missing-braces -Wno-unused-parameter -Wno-unused-variable

$(PSR_TARGETS): %: %.c
	gcc $^ ../src/psr.c $(WARN) -I../inc/ -o $@

$(PSR_AND_PWA_TARGETS): %: %.c
	gcc $^ ../src/psr.c ../src/pwa.c $(WARN) -I../inc/ -lgdi32 -o $@

all: $(PSR_TARGETS) $(PSR_AND_PWA_TARGETS)

.PHONY: clean
clean:
	del *.exe
	del *.bmp
