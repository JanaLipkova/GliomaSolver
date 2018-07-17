/*
 *  MRAGHuffmanEncoder.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 1/9/09.
 *  Copyright 2009 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include "MRAGEncoder.h"
#include <map>
#undef min
#undef max
#include <algorithm>
#include <queue>
#include <stack>


using namespace std;

namespace MRAG
{
	template<typename DataType>
	class HuffmanEncoder: public Encoder<DataType>
	{
		static const bool bDebug = false;
		static const bool bCheat = false;
		static const bool bQComprWhenCheating = true;
		vector<DataType> m_vCheat;
		
		bool m_bPathologicCase;
		DataType* m_PathologicSymbol;
		double m_entropyPerDigit;
		int m_maxDigitEncodingLength;
		
		Encoder<DataType> m_encodedTable_Symbols;
		vector<BitStream> m_encodedTable_Encodings;
		BitStream m_encodedStream;
		
		template<typename Stream>
		void _computeFrequencies(const Stream& stream, map<DataType, float>& mapFrequencies, double& approxEntropy) const
		{
			//1. fill the freq map
			//2. normalize it, computing the probability approximations
			
			//1.
			{
				const typename Stream::const_iterator itEnd = stream.end();
				
				for(typename Stream::const_iterator it = stream.begin(); it != itEnd; it++)
					mapFrequencies[*it] += 0.01f;
			}
			
			//2.
			{
				double H = 0;
				const double factor = 1./(stream.size()*0.01);
				
				const typename map<DataType, float>::iterator itEnd = mapFrequencies.end();
				
				for(typename map<DataType, float>::iterator it = mapFrequencies.begin(); it != itEnd; it++)
				{
					float& freq = it->second;
					
					assert(freq > 0);
					
					freq *= factor;
					H += -freq*log2(freq);
				}
				
				approxEntropy = H;
				
				if(bDebug)
					for(typename map<DataType, float>::iterator it = mapFrequencies.begin(); it != itEnd; it++)
						printf("p[%d] = %f\n", (unsigned int)it->first, it->second);
			}
		}
		
		struct TranslationItem
		{
			float p;
			stack<char> encoding;
			const DataType* symbol;
			
			TranslationItem(): p(0.), encoding(), symbol(){}
			TranslationItem(const TranslationItem& item): p(item.p), encoding(item.encoding), symbol(item.symbol){}
			
		private:
			TranslationItem& operator=(TranslationItem& item){abort(); return *this;}
		};
		
		struct PQitem
		{
			float p;
			vector< TranslationItem* > symbolsAndPrefixes; 
			
			inline bool operator< (const PQitem& item) const
			{
				return (p>item.p);
			}
			
			void add_bit(const int bit) const
			{
				assert(bit==0 || bit==1);
				
				const typename vector< TranslationItem* >::const_iterator itEnd = symbolsAndPrefixes.end();
				for(typename vector< TranslationItem* >::const_iterator it = symbolsAndPrefixes.begin(); it!=itEnd; it++)
					(*it)->encoding.push(bit);
			}
			
			PQitem(float p_,  TranslationItem& item):p(p_),  symbolsAndPrefixes() 
			{
				symbolsAndPrefixes.push_back(&item);
			}
			
			PQitem(const PQitem& item): p(item.p), symbolsAndPrefixes(item.symbolsAndPrefixes){}
			
			PQitem(const PQitem& item1, const PQitem& item2): p(0), symbolsAndPrefixes()
			{
				item1.add_bit(0);
				item2.add_bit(1);
				
				p = item1.p + item2.p;
				
				symbolsAndPrefixes.resize(item1.symbolsAndPrefixes.size()+item2.symbolsAndPrefixes.size());
				set_union(
					item1.symbolsAndPrefixes.begin(), item1.symbolsAndPrefixes.end(),
					item2.symbolsAndPrefixes.begin(), item2.symbolsAndPrefixes.end(),
					symbolsAndPrefixes.begin());
			}
			
		};
		
		void _computeTranslationTable(const map<DataType, float>& mapFrequencies, map<DataType, BitStream >& mapTable, int& maxLength) const
		{
			//1. initialize the translation items
			//2. put them in the priority queue
			//3. merge them, build the tree
			//4. extract the encoding for every symbol
			
			assert(mapTable.size() == 0);
			vector<TranslationItem> translations(mapFrequencies.size());
			
			//1.
			{
				typename map<DataType, float>::const_iterator itSource = mapFrequencies.begin();
				const typename vector<TranslationItem>::iterator itEnd = translations.end();
				for(typename vector<TranslationItem>::iterator it = translations.begin(); it != itEnd; it++, itSource++)
				{
					it->symbol = &itSource->first;
					it->p = itSource->second;
				}
			}
			
			priority_queue<PQitem> mypq;
			
			//2.
			{
				const typename vector<TranslationItem>::iterator itEnd = translations.end();
				for(typename vector<TranslationItem>::iterator it = translations.begin(); it != itEnd; it++)
					mypq.push(PQitem(it->p, *it));
			}
			
			//3.	
			while (mypq.size()>1)
			{
				PQitem i1 = mypq.top();
				mypq.pop();
				PQitem i2 = mypq.top();
				mypq.pop();

				mypq.push(PQitem(i1, i2));
			}
			
			//4.
			{
				int maxL = 0;
				const typename vector<TranslationItem>::iterator itEnd = translations.end();
				for(typename vector<TranslationItem>::iterator it = translations.begin(); it != itEnd; it++)
				{
					TranslationItem& translation = *it;
					stack<char>& encoding = translation.encoding;
					
					BitStream& bits = mapTable[*translation.symbol];
					assert(bits.bits == NULL);
					
					bits.setup(encoding.size());
					maxL = max(maxL, (int)encoding.size());
					
					while(encoding.size()>0)
					{
						bits.append_bits(encoding.top(),1);
						encoding.pop();
					}
				}	
				
				maxLength = maxL;
				
				if(bDebug)
					for(typename vector<TranslationItem>::iterator it = translations.begin(); it != itEnd; it++)
					{
						TranslationItem& translation = *it;
						
						BitStream& bits = mapTable[*translation.symbol];
						printf("Encoding of %d:\n", (unsigned int)(*translation.symbol));
						bits.printBits();
					}	
			}
		}
		
		void _encodeTranslationTable(const map<DataType, BitStream >& mapTable, Encoder<DataType>& encodedSymbols, vector<BitStream>& vEncodings) const
		{
			//0. checks & setup
			//1. fill a vector of symbols, fill the encodings vector
			//2. quantize-encode the vector of symbols
			
			//0.
			assert(vEncodings.size() == 0);
			
			const int nSymbols = mapTable.size();
			vEncodings.resize(nSymbols);
			vector<DataType> vSymbols(nSymbols);
			
			//1.
			{
				vector<BitStream>::iterator itDestEncoding = vEncodings.begin();
				typename vector<DataType>::iterator itDestSymbol = vSymbols.begin();
				const typename map<DataType, BitStream >::const_iterator itEnd = mapTable.end();
				for(typename map<DataType, BitStream >::const_iterator it = mapTable.begin(); it!=itEnd; it++, itDestSymbol++, itDestEncoding++)
				{
					*itDestSymbol = it->first;
					*itDestEncoding = it->second;
				}
			}
			
			//2.
			{
				typename vector<DataType>::const_iterator itSymbol = vSymbols.begin();
				unsigned int maxVal = 0;
				for(int i=0; i<nSymbols; i++,itSymbol++)
					maxVal = max(maxVal, (unsigned int)(*itSymbol));
				
				if (Encoder<DataType>::bVerbose) printf("ENCODING TABLE\n");
				encodedSymbols.encode(vSymbols, maxVal+1);
			}
		}
		
		void _decodeTranslationTable(const Encoder<DataType>& encodedSymbols, const vector<BitStream>& vEncodings,  map<DataType, BitStream >& mapTable) const
		{
			//0. checks, setups
			//1. unpack encodedSymbols into a vector of symbols
			//2. fill the map
			
			//0.
			assert(mapTable.size() == 0);
			const int nSymbols = vEncodings.size();
			vector<DataType> vSymbols;
			
			//1.
			{
				encodedSymbols.decode(vSymbols);
				assert(vSymbols.size() == nSymbols);
			}
			
			//2.
			{
				vector<BitStream>::const_iterator itEncoding = vEncodings.begin();
				typename vector<DataType>::const_iterator itSymbol = vSymbols.begin();
				
				for(int i=0; i<nSymbols; i++, itSymbol++, itEncoding++)
					mapTable[*itSymbol] = *itEncoding;
			}
		}
		
		template<typename Stream>
		void _translate(const double H, const Stream& stream, const map<DataType, BitStream >& mapTable, BitStream& encodedStream) const
		{
			const unsigned int estimated_encoding_bits = (unsigned int)ceil(H*stream.size()*1.25);
			encodedStream.setup(estimated_encoding_bits);
			
			const typename Stream::const_iterator itEnd = stream.end();
			for(typename Stream::const_iterator it = stream.begin(); it != itEnd; it++)
			{
				typename map<DataType, BitStream >::const_iterator itBits = mapTable.find(*it);
				encodedStream.append_bits(itBits->second);
			}
			
			if (bDebug)
				encodedStream.printBits();
		}
		
		template<int chunk_size>
		struct LUT // used in the _translate_back (decompression)
		{
			static const int nEntries = 1<<chunk_size;
			struct Item
			{
				LUT * next_lut;
				const DataType * symbol;
				int encoding_length;
				
				Item(): next_lut(NULL), symbol(NULL), encoding_length(0){}
				Item(const Item& item): next_lut(item.next_lut), symbol(item.symbol), encoding_length(item.encoding_length){}
				
				Item& operator=(const Item item)
				{
					next_lut = item.next_lut;
					symbol = item.symbol;
					encoding_length = item.encoding_length;
					
					return *this;
				}
				
				void printItem() const
				{
					if(symbol != NULL)
						printf("symbol='%c', length = %d\n", *symbol, encoding_length);
					else
					{
						assert(next_lut != NULL);
						printf("symbol=another LUT (0x%x), length = nan\n", next_lut);
					}
				}
			};
			
			Item lut_data[nEntries];
			
			LUT()
			{
				for(int i=0; i<nEntries; i++)
					lut_data[i] = Item();
			}
			
			~LUT()
			{
				for(int i=0; i<nEntries; i++)
				{
					if (lut_data[i].next_lut != NULL)
						delete lut_data[i].next_lut;
					
					lut_data[i].next_lut = NULL;
				}
			}
	
			Item& operator[](unsigned int chunk_data)
			{
				assert(chunk_data<nEntries);
				
				return lut_data[chunk_data];
			}
			
			LUT* access_allocate(unsigned int chunk_data)
			{
				assert(chunk_data<nEntries);
				assert(lut_data[chunk_data].symbol == NULL);
				
				if (lut_data[chunk_data].next_lut == NULL)
					lut_data[chunk_data].next_lut = new LUT();
				
				return lut_data[chunk_data].next_lut;
			}
			
			void printLUT() const 
			{
				printf("LUT (0x%x):\n", this);
				for(int i=0; i<nEntries; i++)
				{
					const unsigned int data = i;
					for(int b=chunk_size-1; b>=0; b--) printf("%d", ((data & (1<<b))>>b));
					printf(": ");
					lut_data[i].printItem();
				}
				
				printf("calling print of next LUTs\n");
				
				for(int i=0; i<nEntries; i++)
					if(lut_data[i].symbol == NULL)
						lut_data[i].next_lut->printLUT();
			}
		};
		
		template<typename Stream>
		void _translate_back(const BitStream& encodedStream, const map<DataType, BitStream >& mapTable, Stream& stream) const
		{
			//1. build the lookup table system (LUT)
			//2. use it for decoding each digit
			
			const unsigned int chunk_size = 8;
			
			//1.
			LUT<chunk_size> startLUT;
			{
				//A. iterate on each (symbol, encoding), building the internal subtrees
				//B. ,, ,, , filling the final leaf
				
				const typename map<DataType, BitStream >::const_iterator itEnd = mapTable.end();
				for(typename map<DataType, BitStream >::const_iterator it = mapTable.begin(); it!=itEnd; it++)
				{
					const DataType& symbol = it->first;
					const BitStream& bits = it->second;
					
					LUT<chunk_size>* currLUT = &startLUT;
					const int nof_chunks = (int)ceil(bits.getLength()/(double)chunk_size);
					
					//A.
					{
						for(int c=0; c<nof_chunks-1; c++)
						{
							const unsigned int chunk_data = bits.get_bits(c*chunk_size, chunk_size);
							currLUT = currLUT->access_allocate(chunk_data);
						}
					}
					
					//B.
					{
						const int remaining_bits = bits.getLength() -  chunk_size*(nof_chunks-1) ;
						const int lshift = chunk_size - remaining_bits;
						const unsigned int prefix = bits.get_bits((unsigned int)(bits.getLength() - remaining_bits), remaining_bits);
						const unsigned int postfixes = 1<<lshift;
						for(unsigned int i=0; i<postfixes; i++)
						{
							typename LUT<chunk_size>::Item& item = (*currLUT)[prefix + (i<<remaining_bits)];
							item.encoding_length = bits.getLength();
							item.symbol = &symbol;
						}
					}
				}
			}
			
			if (bDebug)
				startLUT.printLUT();
			
			//2.
			unsigned int total_encoding_length = encodedStream.getLength();
			
			int start = 0;
			typename Stream::iterator itDest = stream.begin();
			const int n = this->m_nItems;
			for(int i=0; i<n; i++, itDest++)
			{
				int bit_consumed = -1;
				LUT<chunk_size>* currLUT = &startLUT;
				
				int curr_start = start;
				
				do
				{
					assert(curr_start<total_encoding_length);
					
					const unsigned int nof_bits_to_read = min((int)chunk_size, (int)(total_encoding_length - curr_start));
					const unsigned int read_chunk = encodedStream.get_bits(curr_start, nof_bits_to_read);
					
					curr_start += nof_bits_to_read;

					typename LUT<chunk_size>::Item& item = (*currLUT)[read_chunk];
					
					const bool bFound = (item.next_lut == NULL);
					
					if (bFound)
					{
						bit_consumed = item.encoding_length;
						*itDest = *item.symbol;
						break;
					}
					else
						currLUT = item.next_lut;
				}
				while(true);
				
				assert(bit_consumed>0);
				start += bit_consumed;	
			}
		}
		
	public:
		HuffmanEncoder(): Encoder<DataType>(),
			m_vCheat(),
			m_encodedStream(),
			m_bPathologicCase(false),
			m_PathologicSymbol(NULL),
			m_entropyPerDigit(HUGE_VAL),
			m_maxDigitEncodingLength(0),
			m_encodedTable_Symbols(),
			m_encodedTable_Encodings()
		{
		}
		
		~HuffmanEncoder()
		{
			if (m_PathologicSymbol!=NULL)
				delete m_PathologicSymbol;
			
			m_PathologicSymbol = NULL;
			
			const typename vector<BitStream >::iterator itEnd = m_encodedTable_Encodings.end();
			for(typename vector<BitStream >::iterator it = m_encodedTable_Encodings.begin(); it!=itEnd; it++)
				it->dispose();
			
			m_encodedStream.dispose();
		}
		
		void encode(const vector<DataType>& stream, const int nSymbols)
		{
			assert(!this->m_bEncoded);
			
			this->m_nItems = stream.size();
			
			if (bCheat)
			{
				if (bQComprWhenCheating)
					Encoder<DataType>::encode(stream, nSymbols);
				else
				{
					m_vCheat.resize(stream.size());
					copy(stream.begin(), stream.end(), m_vCheat.begin());
				}
				
				this->m_bEncoded = true;
				
				return;
			}
			
			map<DataType, float> mapFrequencies;
			map<DataType, BitStream > mapTranslationTable;
			
			_computeFrequencies(stream, mapFrequencies, m_entropyPerDigit);
			
			m_bPathologicCase = (mapFrequencies.size() == 1);
			
			if (!m_bPathologicCase)
			{
				_computeTranslationTable(mapFrequencies, mapTranslationTable, m_maxDigitEncodingLength);
				_translate(m_entropyPerDigit, stream, mapTranslationTable, m_encodedStream);
				_encodeTranslationTable(mapTranslationTable, m_encodedTable_Symbols, m_encodedTable_Encodings);
				
				const float uncompressedMB = stream.size()*sizeof(DataType)/1024./1024.;
				const double compression_rate = uncompressedMB/getMemorySize();
				
				if (Encoder<DataType>::bVerbose)
					printf("HuffmanEncoder Compression Rate %f  (%.2f KB instead of %.2f KB)\n", compression_rate, getMemorySize()*1024., uncompressedMB*1024.);
				
				/*printf("%s, Compression Rate: %f, est. size %d, actualsize = %d (msg symbols: %d, symbols in alphabet %d)\n",
					   estimated_encoding_bits>=compressed_bits?"Successiful Estimated": "Failed Estimation",
					   compression_rate, compressed_bits/8./1024., uncompressed_bits/8./1024., 
					   estimated_encoding_bits, encodedStream.getLength(), stream.size(), mapTable.size());*/
			}
			else
			{
				m_PathologicSymbol = new DataType();
				*m_PathologicSymbol = mapFrequencies.begin()->first;
			}
			
			this->m_bEncoded = true;
		}
		
		void decode(vector<DataType>& stream) const
		{
			assert(this->m_bEncoded);
			assert(stream.size() == 0);
			
			if (bCheat)
			{
				if (bQComprWhenCheating)
					Encoder<DataType>::decode(stream);
				else
				{
					stream.resize(m_vCheat.size());
					copy(m_vCheat.begin(), m_vCheat.end(), stream.begin());
				}
				
				return;
			}
			
			stream.resize(this->m_nItems);
			
			if (!m_bPathologicCase)
			{
				map<DataType, BitStream > mapTranslationTable;
				_decodeTranslationTable(m_encodedTable_Symbols, m_encodedTable_Encodings, mapTranslationTable);
				_translate_back(m_encodedStream, mapTranslationTable, stream);
			}
			else
			{
				typename vector<DataType>::iterator itDest = stream.begin();
				
				const int n = this->m_nItems;
				for(int i=0; i<n; i++, itDest++)
					*itDest = *m_PathologicSymbol;
			}
		}
		
		float getMemorySize() const
		{
			int sBytes = 0;
			
			if (m_bPathologicCase)
				sBytes+=sizeof(DataType);
			else
			{
				sBytes+= (int)ceil(m_encodedStream.getLength()/8.);
				
				{
					vector<BitStream>::const_iterator itEncoding = m_encodedTable_Encodings.begin();
					const int nEncodings = m_encodedTable_Encodings.size();
					for(int i=0; i<nEncodings; i++,itEncoding++)
						sBytes+= (int)ceil(itEncoding->getLength()/8.);
				}
			}
			
			return (sizeof(HuffmanEncoder) + sBytes)/(1024.*1024.) + m_encodedTable_Symbols.getMemorySize();
		}
		
		
	private:
		//forbidden
		HuffmanEncoder(const HuffmanEncoder& h): Encoder<DataType>(),
			m_encodedStream(),
			m_bPathologicCase(false),
			m_PathologicSymbol(NULL),
		m_entropyPerDigit(HUGE_VAL){ abort(); }
		
		HuffmanEncoder& operator = (const HuffmanEncoder& h){ abort(); return * this;}
	};
}
