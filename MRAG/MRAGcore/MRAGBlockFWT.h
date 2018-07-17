/*
 *  MRAGBlockFWT.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 7/1/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

#include "MRAGHelperMacros.h"
#pragma once
#include "MRAGCommon.h"

#include "MRAGrid.h"
#include "MRAGBlockCollection.h"
#include "MRAGMatrix3D.h"
#include "MRAGBlockLab.h"
#include "MRAGEnvironment.h"
#undef max
#undef min

/*
static bool compute_pos__DEBUG(int _bL, int _bX, int _bY, int &gX, int& gY, bool bChange=false, int iX=-1, int iY=-1)
{
	static int bX = -1;
	static int bY = -1;
	static int bL = -1;
	
	if (bChange)
	{
		bX = _bX;
		bY = _bY;
		bL = _bL;
		return false;
	}
	else
	{
		gX = 8*bX + iX;
		gY = 8*bY + iY;
		
		return (bL==1);
	}
}
*/
namespace MRAG
{

template <int nChannels=1>
class FWTReport
{
private:
	template <typename Wavelets, typename Block , typename Projector,  bool bWriteDetails, int nChannels_,  typename BlockLabType > friend class BlockFWT;
	
	enum FWTReport_Status { FWTReport_Undefined, FWTReport_Open, FWTReport_Concluded};
	
	FWTReport_Status m_status;
	Matrix3D<double>* m_matMinValues[nChannels];
	Matrix3D<double>* m_matMaxValues[nChannels];
	
	FWTReport(int nSizeX, int nSizeY, int nSizeZ):
		m_status(FWTReport_Undefined)
	 { 
		for(int i=0; i<nChannels; i++)
		{
			m_matMinValues[i] = new Matrix3D<double>(nSizeX, nSizeY, nSizeZ);
			m_matMaxValues[i] = new Matrix3D<double>(nSizeX, nSizeY, nSizeZ);
		}
	 }
	
	void clear()
	{
		for(int i=0; i<nChannels; i++)
		{
			*m_matMinValues[i] = (double)HUGE_VAL;
			*m_matMaxValues[i] = 0.0;
		}
	}
	
	void open()
	{
		assert(m_status == FWTReport_Undefined || m_status == FWTReport_Concluded);
		
		clear();
		
		m_status = FWTReport_Open;
	}
	
	template <int iChannel>
	void send(int code, double val)
	{
		assert(m_status == FWTReport_Open);
		
		char idx[3] = {code&1, (code>>1)&1, (code>>2)&1};
		
		double& minVal = m_matMinValues[iChannel]->Access(idx[0], idx[1], idx[2]);
		double& maxVal = m_matMaxValues[iChannel]->Access(idx[0], idx[1], idx[2]);
		
		minVal = std::min(val, minVal);
		maxVal = std::max(val, maxVal);
	}
	
	void conclude()
	{
		assert(m_status == FWTReport_Open);
		
		m_status = FWTReport_Concluded;
	}
	
	const double _getDetailMaxMag(int iChannel) const 
	{
		assert(iChannel>=0 && iChannel<nChannels);
		
		const int n = m_matMaxValues[iChannel]->getNumberOfElements();
		
		double maxVal = m_matMaxValues[iChannel]->LinAccess(1);
		for(int i=2; i<n; i++)
			maxVal = std::max(maxVal, m_matMaxValues[iChannel]->LinAccess(i));
		
		return maxVal;
	}
public:
	FWTReport(const FWTReport& report): 
		m_status(report.m_status)
	{
		const int nSizeX = report.m_matMinValues[0]->getSize()[0];
		const int nSizeY = report.m_matMinValues[0]->getSize()[1];
		const int nSizeZ = report.m_matMinValues[0]->getSize()[2];
		
		for(int i=0; i<nChannels; i++)
		{
			m_matMinValues[i] = new Matrix3D<double>(nSizeX, nSizeY, nSizeZ);
			m_matMaxValues[i] = new Matrix3D<double>(nSizeX, nSizeY, nSizeZ);
			
			*m_matMinValues[i] = *report.m_matMinValues[i];
			*m_matMaxValues[i] = *report.m_matMaxValues[i];
		}
	}
	
	FWTReport& operator=(const FWTReport& report) 
	{
		m_status = report.m_status;
		
		for(int i=0; i<nChannels; i++)
		{
			*m_matMinValues[i] = *report.m_matMinValues[i];
			*m_matMaxValues[i] = *report.m_matMaxValues[i];
		}
		
		return *this;
	}
	
	~FWTReport()
	{
		for(int i=0; i<nChannels; i++)
		{
			delete m_matMinValues[i];
			m_matMinValues[i] = NULL;
			
			delete m_matMaxValues[i]; 
			m_matMaxValues[i] = NULL;
		}
	}
	
	template<int iChannel>
	const double getCoeffMinMag(int ix, int iy=0, int iz=0) const 
	{ 
		assert(m_status == FWTReport_Concluded);
		return m_matMinValues[iChannel]->Read(ix, iy, iz); 
	}
	
	template<int iChannel>
	const double getCoeffMaxMag(int ix, int iy=0, int iz=0) const 
	{ 
		assert(m_status == FWTReport_Concluded);
		return m_matMaxValues[iChannel]->Read(ix, iy, iz); 
	}
	
	template<int iChannel>
	const double getDetailMaxMag() const 
	{
		return _getDetailMaxMag(iChannel);
	}
	
	const double getOverAll_DetailMaxMag() const 
	{
		double maxVal = getDetailMaxMag<0>();
		for(int i=1; i<nChannels; i++)
			maxVal = std::max(maxVal, _getDetailMaxMag(i));
		
		return maxVal;
	}
};

/**
 * Used for automatic refinement/compression.
 * Projector needs to provide a static function Real Project(const T&) which should work for T = Block::ElementType.
 * Easiest way to do this is by using macro #make_projector.
 * @see Science::AutomaticRefinement, Science::AutomaticCompression, MRAGHelperMacros.h, #make_projector
 */
template <typename Wavelets, typename Block , typename Projector = dummy_projector ,bool bWriteDetails=true, int nChannels=1, typename BlockLabType = BlockLab<Block> >
class BlockFWT
{
	typedef Wavelets W;
	typedef typename Block::ElementType Element;
private:
	static const int nReportSizeX = (Block::sizeX>1)?2:1;
	static const int nReportSizeY = (Block::sizeY>1)?2:1;
	static const int nReportSizeZ = (Block::sizeZ>1)?2:1;
	
	enum eBlockFWT_State{eBlockFWT_NotReady, eBlockFWT_Initialized, eBlockFWT_Transformed, eBlockFWT_Reported};
	
	eBlockFWT_State m_state;
	BlockLabType m_blockLab;
	//BlockLab<Block> m_blockLab;
	FWTReport<nChannels> m_report;
	const BlockCollection<Block> * m_refCollection;
	Real * m_weightsHa, * m_weightsGa;
	
	template<bool bWaveletX, bool bWaveletY, bool bWaveletZ, int iChannel>
	Real _filter(const int dx, const int dy, const int dz) const
	{
		const int s[3] = {
			!Block::shouldProcessDirectionY?0 : bWaveletX?Wavelets::GaSupport[0]:Wavelets::HaSupport[0],
			!Block::shouldProcessDirectionY?0 : bWaveletY?Wavelets::GaSupport[0]:Wavelets::HaSupport[0],
			!Block::shouldProcessDirectionZ?0 : bWaveletZ?Wavelets::GaSupport[0]:Wavelets::HaSupport[0]
		};
		
		const int e[3] = {
			!Block::shouldProcessDirectionX?1 : bWaveletX?Wavelets::GaSupport[1]:Wavelets::HaSupport[1],
			!Block::shouldProcessDirectionY?1 : bWaveletY?Wavelets::GaSupport[1]:Wavelets::HaSupport[1],
			!Block::shouldProcessDirectionZ?1 : bWaveletZ?Wavelets::GaSupport[1]:Wavelets::HaSupport[1]
		};
		
		const int o[3] = {
			!Block::shouldProcessDirectionX?0:dx*2,
			!Block::shouldProcessDirectionY?0:dy*2,
			!Block::shouldProcessDirectionZ?0:dz*2
		};
		
		Real result = 0;
		for(int iz=s[2]; iz<e[2]; iz++)
		{
			const Real wZ = (!Block::shouldProcessDirectionZ)?1: (bWaveletZ ? m_weightsGa[iz]:m_weightsHa[iz]);
			for(int iy=s[1]; iy<e[1]; iy++)
			{
				const Real wY = (!Block::shouldProcessDirectionY)?1: (bWaveletY ? m_weightsGa[iy]:m_weightsHa[iy]);
				const Real wZwY = wZ*wY;
				
				const Element * ptr = &m_blockLab.read(o[0] - s[0], o[1] - iy, o[2] - iz);
				
				for(int ix=s[0]; ix<e[0]; ix++, ptr--)
					result += (Projector::template Project<Element, iChannel>(*ptr)) * (wZwY*(bWaveletX?m_weightsGa[ix]:m_weightsHa[ix]));
			}
		}
			
		return result;
	}
	
	template<int iChannel>
	void _fwt(const BlockInfo& info, Block& block)
	{
		/*int a,b;
		compute_pos__DEBUG(info.level, info.index[0], info.index[1], a, b, true);*/
		const int nX = max(1,Block::sizeX/2);
		const int nY = max(1,Block::sizeY/2);
		const int nZ = max(1,Block::sizeZ/2);
		//printf("FST of block %d %d l=%d\n", info.index[0], info.index[1], info.level);
		for(int iz=0; iz<nZ; iz++)
			for(int iy=0; iy<nY; iy++)
				for(int ix=0; ix<nX; ix++)
				{
					Real values[8];
					
					block.wtdata(values[0] = _filter<0,0,0, iChannel>( ix,iy,iz), 0, ix,iy,iz);
					block.wtdata(values[1] = _filter<1,0,0, iChannel>( ix,iy,iz), 1,  ix,iy,iz);
					
					if (Block::shouldProcessDirectionY)
					{
						block.wtdata(values[2] = _filter<0,1,0, iChannel>( ix,iy,iz), 2,  ix,iy,iz);
						block.wtdata(values[3] = _filter<1,1,0, iChannel>( ix,iy,iz), 3,  ix,iy,iz);
						
						if (Block::shouldProcessDirectionZ)
						{
							block.wtdata(values[4] = _filter<0,0,1, iChannel>(ix,iy,iz), 4, ix,iy,iz);
							block.wtdata(values[5] = _filter<1,0,1, iChannel>(ix,iy,iz), 5, ix,iy,iz);
							block.wtdata(values[6] = _filter<0,1,1, iChannel>(ix,iy,iz), 6, ix,iy,iz);
							block.wtdata(values[7] = _filter<1,1,1, iChannel>(ix,iy,iz), 7, ix,iy,iz);
						}
					}
					
					const int n = !Block::shouldProcessDirectionY?2: !Block::shouldProcessDirectionZ? 4 : 8;
					for(int i=0; i<n; i+=2)
					{
						m_report.template send<iChannel>(i, fabs(values[i]));
						m_report.template send<iChannel>(i+1, fabs(values[i+1]));
					}
				}
	}
	
	template <int iFirst, int iLast, int iSize>
	class MultiChannel_FWT
	{
	public:
		inline static void multichannel_fwt(const BlockInfo& info, Block& block,  BlockFWT< Wavelets, Block, Projector, bWriteDetails, nChannels, BlockLabType>& blockfwt)
		{
			blockfwt._fwt<iFirst>(info, block);
			MultiChannel_FWT<iFirst+1, iLast, iSize-1>::multichannel_fwt(info, block, blockfwt);
		}
	};
	
	template <int iFirst, int iLast>
	class MultiChannel_FWT<iFirst, iLast, 0>
	{
	public:
		inline static void multichannel_fwt(const BlockInfo& info, Block& block, BlockFWT< Wavelets, Block, Projector, bWriteDetails, nChannels, BlockLabType>& blockfwt){}
	};
	
public:
	BlockFWT(): 
	m_blockLab(), 
	m_state(eBlockFWT_NotReady),
	m_refCollection(NULL),
	m_report(nReportSizeX, nReportSizeY, nReportSizeZ),
	m_weightsHa(NULL), 
	m_weightsGa(NULL)
	{
		m_weightsGa = new Real[W::GaSupport[1] - W::GaSupport[0]];
		m_weightsHa = new Real[W::HaSupport[1] - W::HaSupport[0]];
		
		//alone in the dark
		
		m_weightsGa -= W::GaSupport[0];
		m_weightsHa -= W::HaSupport[0];
		
		for(int i=W::GaSupport[0]; i<W::GaSupport[1]; i++)
			m_weightsGa[i] = W::getGa(i);
		
		for(int i=W::HaSupport[0]; i<W::HaSupport[1]; i++)
			m_weightsHa[i] = W::getHa(i);
	}
	
	~BlockFWT()
	{
		m_weightsGa += W::GaSupport[0];
		m_weightsHa += W::HaSupport[0];
		
		delete [] m_weightsGa; m_weightsGa = NULL;
		delete [] m_weightsHa; m_weightsHa = NULL;
	}
	
	void prepare(const BlockCollection<Block>& collection, BoundaryInfo& boundaryInfo)
	{
		const int stencilStart[3] = {
			!Block::shouldProcessDirectionX?0: 1+min(-Wavelets::GaSupport[1], -Wavelets::HaSupport[1]),
			!Block::shouldProcessDirectionY?0: 1+min(-Wavelets::GaSupport[1], -Wavelets::HaSupport[1]),
			!Block::shouldProcessDirectionZ?0: 1+min(-Wavelets::GaSupport[1], -Wavelets::HaSupport[1])
		};
		
		const int stencilEnd[3] = {
			!Block::shouldProcessDirectionX?1 : max(-Wavelets::GaSupport[0], -Wavelets::HaSupport[0]),
			!Block::shouldProcessDirectionY?1 : max(-Wavelets::GaSupport[0], -Wavelets::HaSupport[0]),
			!Block::shouldProcessDirectionZ?1 : max(-Wavelets::GaSupport[0], -Wavelets::HaSupport[0])
		};	
		
		m_refCollection = &collection;
		
		m_blockLab.prepare(collection, boundaryInfo, stencilStart, stencilEnd);
		
		m_state = eBlockFWT_Initialized;
	}
	
	const FWTReport<nChannels>& getReport()
	{
		assert(m_state == eBlockFWT_Transformed);
		
		return m_report;
	}
		
	template<int iChannel>
	void fwt(const BlockInfo& info)
	{
		assert(m_state == eBlockFWT_Initialized || m_state == eBlockFWT_Transformed);
		
		Block& block = (*m_refCollection)[info.blockID];
		
		m_report.open();
		m_blockLab.load(info);
		
		_fwt<iChannel>(info, block);
		
		m_state = eBlockFWT_Transformed;
		m_report.conclude();
	}
	
	template<int iFirstChannel, int iLastChannel>
	void multichannel_fwt(const BlockInfo& info)
	{
		assert(m_state == eBlockFWT_Initialized || m_state == eBlockFWT_Transformed);
		
		//Block& block = (*m_refCollection)[info.blockID];
		Block& block = m_refCollection->lock(info.blockID);
		
		m_report.open();
		m_blockLab.load(info);
		
		MultiChannel_FWT<iFirstChannel, iLastChannel, iLastChannel-iFirstChannel+1>::multichannel_fwt(info, block, *this);
		
		m_state = eBlockFWT_Transformed;
		m_report.conclude();

		m_refCollection->release(info.blockID);
	}
	
	template<int iFirstChannel, int iLastChannel>
	struct MultiChannelFWT_Body
	{
		const vector<MRAG::BlockInfo>& vInput;
		vector<FWTReport<nChannels> *> vOutput;
		const BlockCollection<Block>& collection;
		BoundaryInfo& boundaryInfo;
		
		
		MultiChannelFWT_Body(const vector<MRAG::BlockInfo>& vInfo,
							 const BlockCollection<Block>& collection_, BoundaryInfo& boundaryInfo_): 
			vInput(vInfo), vOutput(vInfo.size()), boundaryInfo(boundaryInfo_), collection(collection_)
		{
		}
		
		template<typename BlockedRange>
		void operator()(const BlockedRange& r) const
		{
			BlockFWT fwt;
			
			fwt.prepare(collection, boundaryInfo);

			for(int i=r.begin(); i!=r.end(); i++)
			{
				const BlockInfo& info = vInput[i];
				
				fwt.multichannel_fwt<iFirstChannel, iLastChannel>(info);
				
				*(vOutput[i]) = fwt.getReport();
			}
		}
	};
	
	template<int iFirstChannel, int iLastChannel>
	static vector<FWTReport<nChannels> > multichannel_fwt(
		const vector<MRAG::BlockInfo>& vInfos,
		const BlockCollection<Block>& collection, BoundaryInfo& boundaryInfo)
	{
		const int n = vInfos.size();
		
		MultiChannelFWT_Body<iFirstChannel, iLastChannel> body(vInfos, collection, boundaryInfo);
		
		for(int i=0; i<n; i++)
			body.vOutput[i] = new FWTReport<nChannels>(nReportSizeX, nReportSizeY, nReportSizeZ);
		
#ifdef _MRAG_TBB
		//printf("(BlockFWT::multichannel_fwt: TBB)\n");
		const int nThreads = _MRAG_TBB_NTHREADS_HINT;
		tbb::parallel_for(tbb::blocked_range<size_t>(0, n,std::max(1,n/nThreads)), body);
#else
		body(SimpleInterval(0, n));
#endif
		vector<FWTReport<nChannels> > vResult;
		
		for(int i=0; i<n; i++)
		{
			vResult.push_back(*body.vOutput[i]);
			
			delete body.vOutput[i];
			
			body.vOutput[i] = NULL;
		}
		
		return vResult;
	}
	
private:
	
	//forbidden
	BlockFWT(const BlockFWT&):m_blockLab(), 
	m_state(eBlockFWT_NotReady),
	m_refCollection(NULL),
	m_report(nReportSizeX, nReportSizeY, nReportSizeZ),
	m_weightsHa(NULL), 
	m_weightsGa(NULL){abort();}
	
	BlockFWT& operator=(BlockFWT&){abort(); return *this;}
};

}
