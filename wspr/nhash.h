#pragma once

#include "types.h"

#include <stddef.h>
#include <stdint.h>	/* defines uint32_t etc */

#ifdef __cplusplus
extern "C" {
#endif

u2_t nhash(const void *key, size_t length, uint32_t initval);

#ifdef __cplusplus
}
#endif
