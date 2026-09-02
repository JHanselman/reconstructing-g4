# External C solver for reconstructing-g4:
#   acb_theta_cli -- theta functions and their z-derivatives via FLINT's
#                    acb_theta module (used by ThetaFlint in magma/FlintWrapper.m)
#
# Needs FLINT >= 3.1 (3.2 or later preferred; acb_theta was added in 3.1).
#
# Linux:  make            (needs libflint-dev >= 3.1, libmpfr-dev, libgmp-dev;
#                          note Ubuntu 24.04 ships FLINT 3.0.1, which is too old --
#                          use `make server` there)
# macOS:  make mac        (brew install flint mpfr gmp)
# Test:   make test
#
# The Magma glue looks for the binary at $THETA_FLINT_BIN, then at
# <repo>/Cext/acb_theta_cli, then on $PATH.

CC = gcc
CFLAGS = -O2
LIBS = -lflint -lmpfr -lgmp -lm

FLINT_BREW_FLAGS = \
	-I$$(brew --prefix flint)/include -L$$(brew --prefix flint)/lib \
	-I$$(brew --prefix mpfr)/include -L$$(brew --prefix mpfr)/lib \
	-I$$(brew --prefix gmp)/include  -L$$(brew --prefix gmp)/lib

all: acb_theta_cli

acb_theta_cli: acb_theta_cli.c
	$(CC) $(CFLAGS) acb_theta_cli.c -o acb_theta_cli $(LIBS)

mac: acb_theta_cli.c
	clang -O2 acb_theta_cli.c -o acb_theta_cli $(FLINT_BREW_FLAGS) $(LIBS)

# For servers without admin access or with an old system FLINT: first run
#   sh build_deps.sh
# to build GMP/MPFR/FLINT into PREFIX, then this target static-links against
# them (no LD_LIBRARY_PATH needed at runtime).
PREFIX ?= $(HOME)/.local/acb_theta
server: acb_theta_cli.c
	$(CC) $(CFLAGS) acb_theta_cli.c -o acb_theta_cli -I$(PREFIX)/include \
	  $(PREFIX)/lib/libflint.a $(PREFIX)/lib/libmpfr.a $(PREFIX)/lib/libgmp.a -lm

test: acb_theta_cli
	./acb_theta_cli --selftest

clean:
	rm -f acb_theta_cli

.PHONY: all mac test clean server
