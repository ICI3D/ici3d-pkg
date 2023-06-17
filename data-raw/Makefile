
default: all

R = $(strip Rscript $^ $(1) $@)

VPATH = $(dir $(shell ls -d */))
SRCS := $(notdir $(shell ls **/*.R))
TARGETS := $(addprefix data/,$(SRCS:R=rda))

data/%.rda: %.R
	$(call R)

all: ${TARGETS}
