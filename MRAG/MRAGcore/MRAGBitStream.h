/*
 *  MRAGBitStream.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 1/10/09.
 *  Copyright 2009 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include <cstring>

namespace MRAG
{
	struct BitStream
	{
		unsigned int * bits;
		unsigned int word32_index;
		unsigned int words32;
		unsigned int bits_end;
		
		BitStream(): bits(NULL), words32(0), word32_index(0), bits_end(0){}
		BitStream(const BitStream& b): bits(), words32(b.words32), word32_index(b.word32_index), bits_end(b.bits_end)
		{
			assert (bits == NULL);
			bits = b.bits;
		}
		
		BitStream& operator=(const BitStream& b)
		{
			assert(bits == NULL);
			memcpy(this, &b, sizeof(BitStream));
			return *this;
		}
		
		~BitStream(){}
		
		void dispose()
		{
			if (bits != NULL ) free(bits);
			bits = NULL;
		}
		
		void setup(const int nSizeInBits=0)
		{
			assert(bits == NULL);
			
			const unsigned int nMinimumWords32 = 10;
			words32 = max(nMinimumWords32, (unsigned int)ceil(nSizeInBits/32.));
			bits = (unsigned int *) malloc(sizeof(unsigned int)*words32);
		}
		
		inline void append_bits(unsigned int new_bits, unsigned int nof_bits)
		{
			assert(nof_bits>0);
			assert(bits != NULL);
			
			if (word32_index + (unsigned int)ceil((nof_bits+bits_end-1)/32.) >= words32)
			{
				words32 *= 2;
				bits =  (unsigned int *) realloc(bits, words32*sizeof(unsigned int));
				assert(bits != NULL);
			}

			_CopyToTheEnd(bits[word32_index], new_bits, bits_end);
			
			const int copied_bits = min(32-bits_end, nof_bits);
			
			word32_index += (bits_end + copied_bits & 0x00000020)>>5;
			bits_end = (bits_end + copied_bits) & 0x0000001F;
			
			if (copied_bits < nof_bits)
				_CopyFromTheBeginning(bits[word32_index], new_bits, copied_bits);
			
			bits_end += nof_bits-copied_bits;
		}
		
		inline void append_bits(const BitStream& stream)
		{
			for(int w=0; w<(int)(stream.word32_index)-1; w++)
				append_bits(stream.bits[w], 32);
			
			append_bits(stream.bits[stream.word32_index], stream.bits_end);
		}
		
		inline unsigned int get_bits(unsigned int bit_start, unsigned int nof_bits) const
		{
			assert (nof_bits <= 32);
			assert(bit_start + nof_bits <= word32_index*32 + bits_end);
			
			const int word_index = (int)(bit_start/32.);
			const int bit_index = bit_start % 32;

			const unsigned int s1 = bit_index; 
			const unsigned int e1 = min((unsigned int)32, bit_index+nof_bits); 
			const unsigned int L1 = e1 - s1;
			const unsigned int m1 = ((0xFFFFFFFF << 32-e1) >> 32-L1) << s1;
			const unsigned int w1 = (bits[word_index] & m1) >> s1;
			
			const unsigned int e2 = nof_bits - L1;
			const unsigned int t = 32-e2;
			const unsigned int m2 = ((0xFFFFFFFF << t) >> t);
			const unsigned int w2 = (e2?(bits[word_index+1] & m2):0) << L1;
			
			return w1 | w2;
		}
		
		void printBits()
		{
			printf("printBits (from MSB to LSB, i.e. from last to first):\n");
			printf("#bits = %d\n", word32_index*32+bits_end);
			
			for(int w=word32_index; w>=0; w--)
				for(int b=-1+((w==word32_index)?bits_end:32); b>=0; b--)
					printf("%d", ((bits[w] & (1<<b))>>b));
			
			printf("\n");
		}
		
		const unsigned int getLength() const
		{
			return 32*(word32_index) + bits_end;
		}
		
	private:
		inline void _CopyToTheEnd(unsigned int& ob, unsigned int nb, unsigned int be)
		{
			ob &= (be? (0xFFFFFFFF>>32-be):0x00000000);
			ob |= nb<<be;
		}
		
		inline void _CopyFromTheBeginning(unsigned int& ob, unsigned int nb, int ns)
		{
			ob = nb>>ns;
		}
	};
}