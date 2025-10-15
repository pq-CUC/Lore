# ==============================================================================
# ===================      General Configuration      ==========================
# ==============================================================================

# Compiler settings
CC = gcc
LDFLAGS =

# Base compilation flags
CFLAGS_BASE = -O3 -Wall -Wconversion -Wextra -Wpedantic -Wmissing-prototypes -Wredundant-decls -Wshadow -Wpointer-arith -fomit-frame-pointer -I.
# Define separate compilation flags for each security level
CFLAGS_1 = $(CFLAGS_BASE) -DLORE_LEVEL=1
CFLAGS_2 = $(CFLAGS_BASE) -DLORE_LEVEL=2
CFLAGS_3 = $(CFLAGS_BASE) -DLORE_LEVEL=3

# Target directory
OBJDIR = obj

# VPATH tells make where to find source files
VPATH = .:test:kat

# ==============================================================================
# ===================      Source File Lists      ==============================
# ==============================================================================

# Core library source files
SOURCES_BASE = kem.c pke.c indcpa.c polyvec.c poly.c ntt.c sampler.c reduce.c \
               symmetric-shake.c fips202.c toomcook.c verify.c

# Random number source files (distinguished by purpose)
SOURCES_RANDOM_REAL = randombytes.c
SOURCES_RANDOM_KAT = rng.c

# Test program source files
SOURCES_TEST_MAIN = test_lore.c
SOURCES_SPEED_MAIN = test_speed.c cpucycles.c speed_print.c
SOURCES_KAT_MAIN = PQCgenKAT_pke.c

# ==============================================================================
# ===================      Object File Lists (per level)      ====================
# ==============================================================================

# Generate separate object file lists for each level
OBJECTS_BASE_1 = $(addsuffix _1, $(addprefix $(OBJDIR)/, $(SOURCES_BASE:.c=.o)))
OBJECTS_BASE_2 = $(addsuffix _2, $(addprefix $(OBJDIR)/, $(SOURCES_BASE:.c=.o)))
OBJECTS_BASE_3 = $(addsuffix _3, $(addprefix $(OBJDIR)/, $(SOURCES_BASE:.c=.o)))

OBJECTS_RANDOM_REAL_1 = $(addsuffix _1, $(addprefix $(OBJDIR)/, $(SOURCES_RANDOM_REAL:.c=.o)))
OBJECTS_RANDOM_REAL_2 = $(addsuffix _2, $(addprefix $(OBJDIR)/, $(SOURCES_RANDOM_REAL:.c=.o)))
OBJECTS_RANDOM_REAL_3 = $(addsuffix _3, $(addprefix $(OBJDIR)/, $(SOURCES_RANDOM_REAL:.c=.o)))

OBJECTS_RANDOM_KAT_1 = $(addsuffix _1, $(addprefix $(OBJDIR)/, $(SOURCES_RANDOM_KAT:.c=.o)))
OBJECTS_RANDOM_KAT_2 = $(addsuffix _2, $(addprefix $(OBJDIR)/, $(SOURCES_RANDOM_KAT:.c=.o)))
OBJECTS_RANDOM_KAT_3 = $(addsuffix _3, $(addprefix $(OBJDIR)/, $(SOURCES_RANDOM_KAT:.c=.o)))

OBJECTS_TEST_MAIN_1 = $(addsuffix _1, $(addprefix $(OBJDIR)/, $(SOURCES_TEST_MAIN:.c=.o)))
OBJECTS_TEST_MAIN_2 = $(addsuffix _2, $(addprefix $(OBJDIR)/, $(SOURCES_TEST_MAIN:.c=.o)))
OBJECTS_TEST_MAIN_3 = $(addsuffix _3, $(addprefix $(OBJDIR)/, $(SOURCES_TEST_MAIN:.c=.o)))

OBJECTS_SPEED_MAIN_1 = $(addsuffix _1, $(addprefix $(OBJDIR)/, $(SOURCES_SPEED_MAIN:.c=.o)))
OBJECTS_SPEED_MAIN_2 = $(addsuffix _2, $(addprefix $(OBJDIR)/, $(SOURCES_SPEED_MAIN:.c=.o)))
OBJECTS_SPEED_MAIN_3 = $(addsuffix _3, $(addprefix $(OBJDIR)/, $(SOURCES_SPEED_MAIN:.c=.o)))

OBJECTS_KAT_MAIN_1 = $(addsuffix _1, $(addprefix $(OBJDIR)/, $(SOURCES_KAT_MAIN:.c=.o)))
OBJECTS_KAT_MAIN_2 = $(addsuffix _2, $(addprefix $(OBJDIR)/, $(SOURCES_KAT_MAIN:.c=.o)))
OBJECTS_KAT_MAIN_3 = $(addsuffix _3, $(addprefix $(OBJDIR)/, $(SOURCES_KAT_MAIN:.c=.o)))


# ==============================================================================
# ===================      Build Rules (Recipes)      ===========================
# ==============================================================================

# Default target: compile all programs for all levels
all: test speed kat

# Define convenience targets for building different types of programs
test: test_lore_1 test_lore_2 test_lore_3
speed: test_speed_1 test_speed_2 test_speed_3
kat: PQCgenKAT_pke_1 PQCgenKAT_pke_2 PQCgenKAT_pke_3

# --- Linking Rules (per level) ---

# Test programs use the REAL randombytes source
test_lore_1: $(OBJECTS_BASE_1) $(OBJECTS_RANDOM_REAL_1) $(OBJECTS_TEST_MAIN_1)
	@echo "==> Linking $@ ..." && $(CC) $(CFLAGS_1) $^ -o $@ $(LDFLAGS)
test_lore_2: $(OBJECTS_BASE_2) $(OBJECTS_RANDOM_REAL_2) $(OBJECTS_TEST_MAIN_2)
	@echo "==> Linking $@ ..." && $(CC) $(CFLAGS_2) $^ -o $@ $(LDFLAGS)
test_lore_3: $(OBJECTS_BASE_3) $(OBJECTS_RANDOM_REAL_3) $(OBJECTS_TEST_MAIN_3)
	@echo "==> Linking $@ ..." && $(CC) $(CFLAGS_3) $^ -o $@ $(LDFLAGS)

# Speed tests also use the REAL randombytes source
test_speed_1: $(OBJECTS_BASE_1) $(OBJECTS_RANDOM_REAL_1) $(OBJECTS_SPEED_MAIN_1)
	@echo "==> Linking $@ ..." && $(CC) $(CFLAGS_1) $^ -o $@ $(LDFLAGS)
test_speed_2: $(OBJECTS_BASE_2) $(OBJECTS_RANDOM_REAL_2) $(OBJECTS_SPEED_MAIN_2)
	@echo "==> Linking $@ ..." && $(CC) $(CFLAGS_2) $^ -o $@ $(LDFLAGS)
test_speed_3: $(OBJECTS_BASE_3) $(OBJECTS_RANDOM_REAL_3) $(OBJECTS_SPEED_MAIN_3)
	@echo "==> Linking $@ ..." && $(CC) $(CFLAGS_3) $^ -o $@ $(LDFLAGS)

# KAT generation uses the DETERMINISTIC rng source
PQCgenKAT_pke_1: $(OBJECTS_BASE_1) $(OBJECTS_RANDOM_KAT_1) $(OBJECTS_KAT_MAIN_1)
	@echo "==> Linking $@ ..." && $(CC) $(CFLAGS_1) $^ -o $@ $(LDFLAGS)
PQCgenKAT_pke_2: $(OBJECTS_BASE_2) $(OBJECTS_RANDOM_KAT_2) $(OBJECTS_KAT_MAIN_2)
	@echo "==> Linking $@ ..." && $(CC) $(CFLAGS_2) $^ -o $@ $(LDFLAGS)
PQCgenKAT_pke_3: $(OBJECTS_BASE_3) $(OBJECTS_RANDOM_KAT_3) $(OBJECTS_KAT_MAIN_3)
	@echo "==> Linking $@ ..." && $(CC) $(CFLAGS_3) $^ -o $@ $(LDFLAGS)

# --- Compilation Rules (per level) ---

# Rule to create the obj directory
$(OBJDIR):
	@mkdir -p $(OBJDIR)

# Compile Level 1 .o files
$(OBJDIR)/%.o_1: %.c | $(OBJDIR)
	@echo "Compiling $< (Level 1)..."
	@$(CC) $(CFLAGS_1) -c $< -o $@

# Compile Level 2 .o files
$(OBJDIR)/%.o_2: %.c | $(OBJDIR)
	@echo "Compiling $< (Level 2)..."
	@$(CC) $(CFLAGS_2) -c $< -o $@

# Compile Level 3 .o files
$(OBJDIR)/%.o_3: %.c | $(OBJDIR)
	@echo "Compiling $< (Level 3)..."
	@$(CC) $(CFLAGS_3) -c $< -o $@

# --- Clean Rule ---
clean:
	@echo "Cleaning..."
	@-rm -rf $(OBJDIR)
	@-rm -f test_lore_* test_speed_* PQCgenKAT_pke_*
	@-rm -f kat/*.req kat/*.rsp
	@echo "Done."

# Declare these targets as "phony" as they do not represent actual files
.PHONY: all test speed kat clean