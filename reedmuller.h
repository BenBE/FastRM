#pragma once

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

bool reedmuller_isvalid(uint8_t r, uint8_t m);

size_t reedmuller_rawsize(uint8_t r, uint8_t m);
size_t reedmuller_blocksize(uint8_t r, uint8_t m);
size_t reedmuller_maxerror(uint8_t r, uint8_t m);

bool reedmuller_encode(uint8_t r, uint8_t m, const uint8_t* ibuf, size_t isize, uint8_t* obuf, size_t osize);
bool reedmuller_decode(uint8_t r, uint8_t m, uint8_t* ibuf, size_t isize, uint8_t* obuf, size_t osize);
