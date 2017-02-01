/*
    damapper_bwt
    Copyright (C) 2016 German Tischler

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
//#define LASTOBAM_PRINT_SAM

#include <config.h>
#include <libmaus2/aio/SynchronousGenericInput.hpp>
#include <libmaus2/bambam/BamAlignmentEncoderBase.hpp>
#include <libmaus2/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus2/bambam/BamHeader.hpp>
#include <libmaus2/bambam/BamAlignmentDecoder.hpp>
#include <libmaus2/bambam/CigarStringParser.hpp>
#include <libmaus2/bitio/CompactArrayWriterFile.hpp>
#include <libmaus2/dazzler/align/AlignmentWriter.hpp>
#include <libmaus2/dazzler/align/LASToBamConverterBase.hpp>
#include <libmaus2/dazzler/align/OverlapIndexer.hpp>
#include <libmaus2/dazzler/align/SimpleOverlapParser.hpp>
#include <libmaus2/dazzler/align/SortingOverlapOutputBuffer.hpp>
#include <libmaus2/fastx/CoordinateCache.hpp>
#include <libmaus2/fastx/FastaBPDecoder.hpp>
#include <libmaus2/fastx/FastABPGenerator.hpp>
#include <libmaus2/fastx/FastAIndex.hpp>
#include <libmaus2/fastx/LineBufferFastAReader.hpp>
#include <libmaus2/fm/SampledSA.hpp>
#include <libmaus2/lcs/AlignerFactory.hpp>
#include <libmaus2/lcs/DalignerLocalAlignment.hpp>
#include <libmaus2/lcs/LocalAlignmentPrint.hpp>
#include <libmaus2/lcs/NP.hpp>
#include <libmaus2/lz/BgzfDeflate.hpp>
#include <libmaus2/math/numbits.hpp>
#include <libmaus2/parallel/NumCpus.hpp>
#include <libmaus2/rank/DNARankKmerCache.hpp>
#include <libmaus2/timing/RealTimeClock.hpp>
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/util/MemUsage.hpp>
#include <libmaus2/util/OutputFileNameTools.hpp>
#include <libmaus2/util/PrefixSums.hpp>
#include <libmaus2/clustering/HashCreationBase.hpp>
#include <libmaus2/lcs/SimdX86BandedGlobalAlignmentScoreY256_16.hpp>
#include <libmaus2/lcs/NNP.hpp>
#include <libmaus2/fastx/FastAToCompact4BigBand.hpp>
#include <libmaus2/suffixsort/bwtb3m/BwtMergeSort.hpp>
#include <libmaus2/util/FiniteSizeHeap.hpp>
#include <libmaus2/geometry/RangeSet.hpp>
#include <libmaus2/util/GrowingFreeList.hpp>
#include <libmaus2/sorting/InterleavedRadixSort.hpp>
#include <libmaus2/geometry/RangeSet.hpp>
#include <libmaus2/dazzler/align/SimpleOverlapVectorParser.hpp>
#include <libmaus2/math/numbits.hpp>
#include <libmaus2/dazzler/align/RefMapEntryVector.hpp>

#include <libmaus2/lcs/EnvelopeFragment.hpp>
#include <libmaus2/lcs/FragmentEnvelope.hpp>
#include <libmaus2/geometry/RangeSet.hpp>
#include <libmaus2/math/binom.hpp>
#include <libmaus2/bitio/CompactArrayWriter.hpp>
#include <libmaus2/util/WriteableString.hpp>
#include <libmaus2/fastx/StreamFastAReader.hpp>
#include <libmaus2/util/MemoryStatistics.hpp>
#include <libmaus2/fastx/FastAIndexGenerator.hpp>
#include <libmaus2/parallel/PosixSemaphore.hpp>
#include <libmaus2/parallel/SimpleThreadWorkPackage.hpp>
#include <libmaus2/parallel/SimpleThreadWorkPackageDispatcher.hpp>
#include <libmaus2/parallel/SimpleThreadPool.hpp>
#include <libmaus2/util/Demangle.hpp>
#include <libmaus2/bambam/BamAlignmentDecoderFactory.hpp>
#include <libmaus2/bambam/ProgramHeaderLineSet.hpp>
#include <libmaus2/lz/BgzfDeflateOutputBufferBase.hpp>
#include <libmaus2/parallel/LockedGrowingFreeList.hpp>

static int getDefaultLevel() { return Z_DEFAULT_COMPRESSION; }

struct BgzfDeflateOutputBufferBaseTypeInfo
{
	typedef libmaus2::lz::BgzfDeflateOutputBufferBase element_type;
	typedef element_type::shared_ptr_type pointer_type;

	static pointer_type getNullPointer()
	{
		pointer_type p;
		return p;
	}

	static pointer_type deallocate(pointer_type /* p */)
	{
		return getNullPointer();
	}
};

struct BgzfDeflateOutputBufferBaseAllocator
{
	typedef libmaus2::lz::BgzfDeflateOutputBufferBase element_type;

	int level;

	BgzfDeflateOutputBufferBaseAllocator(int const rlevel = Z_DEFAULT_COMPRESSION) : level(rlevel)
	{

	}

	element_type::shared_ptr_type operator()() const
	{
		element_type::shared_ptr_type tptr(new element_type(level));
		return tptr;
	}
};

struct BgzfDeflateZStreamBaseTypeInfo
{
	typedef libmaus2::lz::BgzfDeflateZStreamBase element_type;
	typedef element_type::shared_ptr_type pointer_type;

	static pointer_type getNullPointer()
	{
		pointer_type p;
		return p;
	}

	static pointer_type deallocate(pointer_type /* p */)
	{
		return getNullPointer();
	}
};

struct BgzfDeflateZStreamBaseAllocator
{
	typedef libmaus2::lz::BgzfDeflateZStreamBase element_type;

	int level;

	BgzfDeflateZStreamBaseAllocator(int const rlevel = Z_DEFAULT_COMPRESSION) : level(rlevel)
	{

	}

	element_type::shared_ptr_type operator()() const
	{
		element_type::shared_ptr_type tptr(new element_type(level));
		return tptr;
	}
};

// damapper interface

#if defined(__cplusplus)
extern "C" {
#endif

#include "QV.h"
#include "DB.h"
#include "align.h"
#include "map.h"

#if defined(__cplusplus)
}
#endif

uint64 MEM_LIMIT = 0;
uint64 MEM_PHYSICAL = 0;
int VERBOSE = 0;
int BIASED = 0;
int SPACING = 100;
int PROFILE = 0;
double BEST_TIE = 1.0;

struct AlignerAllocator
{
	AlignerAllocator() {}

	libmaus2::lcs::Aligner::shared_ptr_type operator()() const
	{
		std::set<libmaus2::lcs::AlignerFactory::aligner_type> const S = libmaus2::lcs::AlignerFactory::getSupportedAligners();

		if ( S.find(libmaus2::lcs::AlignerFactory::libmaus2_lcs_AlignerFactory_y256_8) != S.end() )
		{
			libmaus2::lcs::Aligner::unique_ptr_type T(libmaus2::lcs::AlignerFactory::construct(libmaus2::lcs::AlignerFactory::libmaus2_lcs_AlignerFactory_y256_8));
			libmaus2::lcs::Aligner::shared_ptr_type S(T.release());
			return S;
		}
		else if ( S.find(libmaus2::lcs::AlignerFactory::libmaus2_lcs_AlignerFactory_x128_8) != S.end() )
		{
			libmaus2::lcs::Aligner::unique_ptr_type T(libmaus2::lcs::AlignerFactory::construct(libmaus2::lcs::AlignerFactory::libmaus2_lcs_AlignerFactory_x128_8));
			libmaus2::lcs::Aligner::shared_ptr_type S(T.release());
			return S;
		}
		else if ( S.find(libmaus2::lcs::AlignerFactory::libmaus2_lcs_AlignerFactory_NP) != S.end() )
		{
			libmaus2::lcs::Aligner::unique_ptr_type T(libmaus2::lcs::AlignerFactory::construct(libmaus2::lcs::AlignerFactory::libmaus2_lcs_AlignerFactory_x128_8));
			libmaus2::lcs::Aligner::shared_ptr_type S(T.release());
			return S;
		}
		else if ( S.size() )
		{
			libmaus2::lcs::Aligner::unique_ptr_type T(libmaus2::lcs::AlignerFactory::construct(*(S.begin())));
			libmaus2::lcs::Aligner::shared_ptr_type S(T.release());
			return S;
		}
		else
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "AlignerAllocator: no aligners found" << std::endl;
			lme.finish();
			throw lme;
		}

	}
};

struct OverlapComparatorBReadRefLength
{
	bool operator()(libmaus2::dazzler::align::Overlap const & lhs, libmaus2::dazzler::align::Overlap const & rhs) const
	{
		if ( lhs.bread != rhs.bread )
			return lhs.bread < rhs.bread;

		return (lhs.path.aepos-lhs.path.abpos) > (rhs.path.aepos-rhs.path.abpos);
	}
};

#if 0
void writeDatabaseAsFasta(HITS_DB const & db, std::ostream & out, std::vector<std::string> const * readnames = 0)
{
	uint64_t const towrite =
		std::min(
			static_cast<uint64_t>(db.nreads),
			static_cast<uint64_t>(readnames ? readnames->size() : db.nreads)
		);
	for ( uint64_t i = 0; i < towrite; ++i )
	{
		std::ostringstream namestr;

		if ( readnames )
			namestr << ((*readnames)[i]);
		else
			namestr << "L0/" << i << "/" << 0 << "_" << db.reads[i].rlen << " RQ=0.851";

		out << ">" << namestr.str() << '\n';
		char const * p = reinterpret_cast<char const *>(db.bases) + db.reads[i].boff;
		uint64_t n = db.reads[i].rlen;
		uint64_t const colwidth = 80;

		while ( n )
		{
			uint64_t const towrite = std::min(colwidth,n);
			for ( uint64_t i = 0; i < towrite; ++i )
				out.put(libmaus2::fastx::remapChar(p[i]));
			out.put('\n');
			p += towrite;
			n -= towrite;
		}
	}
}

void writeDatabaseAsFasta(HITS_DB const & db, std::string const & fn, std::vector<std::string> const & readnames)
{
	libmaus2::aio::OutputStreamInstance::unique_ptr_type OSI(new libmaus2::aio::OutputStreamInstance(fn+".fasta"));
	writeDatabaseAsFasta(db,*OSI,NULL);
	OSI->flush();
	OSI.reset();

	std::ostringstream delstr;
	delstr << "/home/tischler/src/git/DAZZ_DB/DBrm " << fn << ".db";
	std::string const del = delstr.str();
	system(del.c_str());

	std::ostringstream comstr;
	comstr << "/home/tischler/src/git/DAZZ_DB/fasta2DB " << fn << ".db " << fn << ".fasta";
	std::string const com = comstr.str();
	system(com.c_str());

	std::ostringstream splitcomstr;
	splitcomstr << "/home/tischler/src/git/DAZZ_DB/DBsplit -s500 -x0 " << fn << ".db";
	std::string const splitcom = splitcomstr.str();
	system(splitcom.c_str());

	libmaus2::aio::OutputStreamInstance::unique_ptr_type OSIrn(new libmaus2::aio::OutputStreamInstance(fn+".fasta"));
	writeDatabaseAsFasta(db,*OSIrn,&readnames);
	OSIrn->flush();
	OSIrn.reset();

	std::ostringstream faidxcomstr;
	faidxcomstr << "samtools faidx " << fn << ".fasta";
	std::string const faidxcom = faidxcomstr.str();
	system(faidxcom.c_str());
}
#endif

struct DalignerInfo
{
	typedef DalignerInfo this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	Work_Data * wd;
	Align_Spec * aspec;
	::Alignment align;
	::Path path;

	void cleanup()
	{
		if ( wd )
		{
			Free_Work_Data(wd);
			wd = 0;
		}
		if ( aspec )
		{
			Free_Align_Spec(aspec);
			aspec = 0;
		}
	}

	DalignerInfo(double const rcor = 0.85, int64_t const rtspace = 100)
	: wd(0), aspec(0)
	{
		float freq[4] = { 0.25, 0.25, 0.25, 0.25 };
		wd = ::New_Work_Data();
		if ( ! wd )
		{
			cleanup();
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "DalignerInfo: failed call New_Work_Data" << std::endl;
			lme.finish();
			throw lme;
		}

		aspec = New_Align_Spec(rcor,rtspace,&freq[0]);

		if ( ! aspec )
		{
			cleanup();
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "DalignerInfo: failed call New_Align_Spec" << std::endl;
			lme.finish();
			throw lme;
		}
	}

	~DalignerInfo()
	{
		cleanup();
	}
};

struct DNAIndexBase
{
	static uint64_t readSamplingRate(std::string const & fn)
	{
		libmaus2::aio::InputStreamInstance ISI(fn);
		return libmaus2::util::NumberSerialisation::deserialiseNumber(ISI);
	}

	static libmaus2::bitio::CompactArray::unique_ptr_type loadCompactArray(std::string const & fn)
	{
		libmaus2::aio::InputStreamInstance ISI(fn);
		libmaus2::bitio::CompactArray::unique_ptr_type CA(new libmaus2::bitio::CompactArray(ISI));
		return UNIQUE_PTR_MOVE(CA);
	}

	static void generateRCFasta(std::string const & fastaname, std::string const & rcfastaname)
	{
		if ( (! libmaus2::util::GetFileSize::fileExists(rcfastaname)) || (libmaus2::util::GetFileSize::isOlder(rcfastaname,fastaname)) )
		{
			libmaus2::aio::InputStreamInstance ISI(fastaname);
			libmaus2::fastx::StreamFastAReaderWrapper SFA(ISI);
			libmaus2::aio::OutputStreamInstance OSI(rcfastaname);
			libmaus2::fastx::StreamFastAReaderWrapper::pattern_type pattern;

			while ( SFA.getNextPatternUnlocked(pattern) )
			{
				std::string & spattern = pattern.spattern;
				std::reverse(spattern.begin(),spattern.end());
				for ( uint64_t i = 0; i < spattern.size(); ++i )
					spattern[i] = libmaus2::fastx::invertUnmapped(spattern[i]);
				OSI << pattern;
			}
		}
	}

	static void generateCompactRCFromCompact(
		std::string const & compactname, std::string const & compactrcname,
		libmaus2::fastx::FastAIndex const * PFAI,
		int verbose
	)
	{
		// compute compact rc file
		if ( (! libmaus2::util::GetFileSize::fileExists(compactrcname)) || (libmaus2::util::GetFileSize::isOlder(compactrcname,compactname)) )
		{
			if ( verbose > 1 )
				std::cerr << "[V] computing " << compactrcname << " from " << compactname << "...";

			libmaus2::bitio::CompactDecoderWrapper CDW(compactname);
			CDW.seekg(0,std::ios::end);
			uint64_t const n = CDW.tellg();
			CDW.clear();
			CDW.seekg(0,std::ios::beg);
			libmaus2::bitio::CompactArrayWriter CW(compactrcname,n,2);
			libmaus2::autoarray::AutoArray<char> C;

			for ( uint64_t i = 0; i < PFAI->size(); ++i )
			{
				uint64_t const l = (*PFAI)[i].length;
				C.ensureSize(l);
				CDW.read(C.begin(),l);
				std::reverse(C.begin(),C.begin()+l);
				for ( uint64_t j = 0; j < l; ++j )
					C[j] ^= 3;
				CW.write(C.begin(),l);
			}

			if ( verbose > 1 )
				std::cerr << "done." << std::endl;
		}

	}

	static void generateDNARank(std::string const & bwtfn, std::string const & dnarankname, uint64_t const numthreads, int const verbose)
	{
		if ( (! libmaus2::util::GetFileSize::fileExists(dnarankname)) || (libmaus2::util::GetFileSize::isOlder(dnarankname,bwtfn) ) )
		{
			if ( verbose > 0 )
				std::cerr << "[V] creating " << dnarankname << "...";

			libmaus2::rank::DNARank::unique_ptr_type Prank(libmaus2::rank::DNARank::loadFromRunLength(bwtfn,numthreads));
			libmaus2::aio::OutputStreamInstance::unique_ptr_type OSI(new libmaus2::aio::OutputStreamInstance(dnarankname));
			Prank->serialise(*OSI);
			OSI->flush();
			OSI.reset();

			if ( verbose > 0 )
				std::cerr << "done." << std::endl;
		}

	}

	static void generateCompactSA(std::string const & saname, std::string const & compactsaname, std::string const & compactsametaname, int const verbose, uint64_t const n, unsigned int const nbits)
	{
		// construct compact sampled suffix array from non compact version if not done already
		if (
			(! libmaus2::util::GetFileSize::fileExists(compactsaname))
			||
			(libmaus2::util::GetFileSize::isOlder(compactsaname,saname))
		)
		{
			if ( verbose )
				std::cerr << "[V] rewriting sampled suffix array...";
			libmaus2::aio::InputStreamInstance::unique_ptr_type SAistr(new libmaus2::aio::InputStreamInstance(saname));
			libmaus2::aio::SynchronousGenericInput<uint64_t> SGI(*SAistr,8*1024);
			uint64_t const sasamplingrate = SGI.get();
			uint64_t const numsamples = (n + sasamplingrate - 1)/sasamplingrate;
			uint64_t const exptsamples = SGI.get();
			assert ( exptsamples == numsamples );
			uint64_t v = 0;
			libmaus2::bitio::CompactArrayWriterFile::unique_ptr_type CAWF(new libmaus2::bitio::CompactArrayWriterFile(compactsaname,nbits));

			for ( uint64_t i = 0; i < numsamples; ++i )
			{
				bool const ok = SGI.getNext(v);
				assert ( ok );
				CAWF->put(v);
			}

			CAWF->flush();
			CAWF.reset();

			libmaus2::aio::OutputStreamInstance::unique_ptr_type OSI(new libmaus2::aio::OutputStreamInstance(compactsametaname));
			libmaus2::util::NumberSerialisation::serialiseNumber(*OSI,sasamplingrate);
			OSI->flush();
			OSI.reset();

			if ( verbose )
				std::cerr << "done, rate " << sasamplingrate << std::endl;
		}

	}

	static uint64_t generateBP(std::string const & replname, std::string const & bpname, int const verbose)
	{
		if ( (! libmaus2::util::GetFileSize::fileExists(bpname)) || (libmaus2::util::GetFileSize::isOlder(bpname,replname) ) )
		{
			if ( verbose > 0 )
				std::cerr << "[V] creating " << bpname << "...";
			libmaus2::aio::InputStreamInstance::unique_ptr_type replISI(new libmaus2::aio::InputStreamInstance(replname));
			libmaus2::aio::OutputStreamInstance::unique_ptr_type bpOSI(new libmaus2::aio::OutputStreamInstance(bpname));
			uint64_t const s = libmaus2::fastx::FastABPGenerator::fastAToFastaBP(*replISI,*bpOSI,&std::cerr,64*1024 /* bs */);
			bpOSI->flush();
			bpOSI.reset();
			replISI.reset();
			if ( verbose > 0 )
				std::cerr << "done." << std::endl;

			return s;
		}
		else
		{
			return libmaus2::util::GetFileSize::getFileSize(bpname);
		}
	}

	static void generateReplAndCompact(std::string const & fastaname, std::string const & replname, std::string const & compactname)
	{
		// generate modified and compact fasta files if they do not exist
		if (
			(! libmaus2::util::GetFileSize::fileExists(replname))
			|| (libmaus2::util::GetFileSize::isOlder(replname,fastaname) )
			|| (libmaus2::util::GetFileSize::isOlder(compactname,fastaname) )
		)
		{
			int const err = libmaus2::fastx::FastAToCompact4BigBand::fastaToCompact4BigBand(
				std::vector<std::string>(1,fastaname),
				false /* rc */,
				false /* replrc */,
				&(std::cerr),
				compactname
			);
			if ( err != EXIT_SUCCESS )
			{
				libmaus2::exception::LibMausException lme;
				lme.getStream() << "[E] Failed to generate compact FastA file" << std::endl;
				lme.finish();
				throw lme;
			}
		}
	}

	static void generateBwt(
		std::string const & compactname,
		std::string const & bwtfn,
		size_t const constrsasamplingrate,
		size_t const iconstrsasamplingrate,
		uint64_t const bwtconstrmem,
		unsigned int const numthreads
	)
	{
		// compute BWT and sampled suffix array if they do not exist
		if ( (! libmaus2::util::GetFileSize::fileExists(bwtfn)) || (libmaus2::util::GetFileSize::isOlder(bwtfn,compactname) ) )
		{
			std::cerr << "[V] sasamplingrate=" << constrsasamplingrate << std::endl;

			libmaus2::suffixsort::bwtb3m::BwtMergeSortOptions options(
				compactname,
				bwtconstrmem,
				numthreads,
				"compactstream",
				true /* bwtonly */,
				bwtfn + "_tmp",
				std::string(), // sparse
				bwtfn,
				iconstrsasamplingrate /* isa */,
				constrsasamplingrate /* sa */
			);

			libmaus2::suffixsort::bwtb3m::BwtMergeSortResult res = libmaus2::suffixsort::bwtb3m::BwtMergeSort::computeBwt(options,&std::cerr);

			res.computeSampledSuffixArray(constrsasamplingrate /* sa */, iconstrsasamplingrate /* isa */,bwtfn + "_tmp",true /* copy input to memory */,
				numthreads,
				bwtconstrmem, /* max sort mem */
				64, /* max tmp files */
				&(std::cerr)
			);
		}
	}

	static uint64_t rewriteFAI(std::ostream & out, std::string const & fainame)
	{
		libmaus2::fastx::FastAIndex::unique_ptr_type tPFAI(libmaus2::fastx::FastAIndex::load(fainame));
		return tPFAI->serialise(out);
	}

	static uint64_t rewriteGeneric(std::ostream & out, std::string const & fn)
	{
		uint64_t const s = libmaus2::util::GetFileSize::getFileSize(fn);
		libmaus2::aio::InputStreamInstance ISI(fn);
		libmaus2::util::GetFileSize::copy(ISI,out,s);
		return s;
	}

	static uint64_t rewriteDNARank(std::ostream & out, std::string const & dnarankname)
	{
		return rewriteGeneric(out,dnarankname);
	}

	static uint64_t rewriteDNAIndexMetaDataBigBand(std::ostream & out, std::string const & index)
	{
		return rewriteGeneric(out,index);
	}

	static uint64_t rewriteCompactSampledSuffixArray(std::ostream & out, std::string const & meta, std::string const & fn)
	{
		uint64_t s = 0;
		uint64_t const sasamplingrate = readSamplingRate(meta);
		s += libmaus2::util::NumberSerialisation::serialiseNumber(out,sasamplingrate);
		s += rewriteGeneric(out,fn);
		return s;
	}

	static uint64_t rewriteBP(std::ostream & out, std::string const & bp)
	{
		uint64_t s = 0;
		uint64_t const fs = libmaus2::util::GetFileSize::getFileSize(bp);
		s += libmaus2::util::NumberSerialisation::serialiseNumber(out,fs);
		s += rewriteGeneric(out,bp);
		return s;
	}
};

struct DNAIndexBuild : public DNAIndexBase
{
	std::string fastaname;
	std::string rcfastaname;
	std::string fainame;
	std::string basefn;
	std::string replname;
	std::string compactname;
	std::string metaname;
	std::string bwtfn;
	std::string bpname;
	std::string saname;
	std::string compactsaname;
	std::string compactsametaname;
	std::string dnarankname;

	std::string histfn;
	std::string isafn;
	std::string preisafn;
	std::string preisametafn;

	void cleanup() const
	{
		libmaus2::aio::FileRemoval::removeFile(rcfastaname);
		libmaus2::aio::FileRemoval::removeFile(fainame);
		libmaus2::aio::FileRemoval::removeFile(replname);
		libmaus2::aio::FileRemoval::removeFile(compactname);
		libmaus2::aio::FileRemoval::removeFile(metaname);
		libmaus2::aio::FileRemoval::removeFile(bwtfn);
		libmaus2::aio::FileRemoval::removeFile(bpname);
		libmaus2::aio::FileRemoval::removeFile(saname);
		libmaus2::aio::FileRemoval::removeFile(compactsaname);
		libmaus2::aio::FileRemoval::removeFile(compactsametaname);
		libmaus2::aio::FileRemoval::removeFile(dnarankname);
		libmaus2::aio::FileRemoval::removeFile(histfn);
		libmaus2::aio::FileRemoval::removeFile(isafn);
		libmaus2::aio::FileRemoval::removeFile(preisafn);
		libmaus2::aio::FileRemoval::removeFile(preisametafn);
	}

	uint64_t n;
	unsigned int nbits;

	libmaus2::fastx::FastAIndex::unique_ptr_type PFAI;

	void loadFAI(std::string const & fainame, int const verbose)
	{
		if ( verbose > 1 )
			std::cerr << "[V] loading FastA (.fai) index...";
		libmaus2::fastx::FastAIndex::unique_ptr_type tPFAI(libmaus2::fastx::FastAIndex::load(fainame));
		PFAI = UNIQUE_PTR_MOVE(tPFAI);
		if ( verbose > 1 )
			std::cerr << "done." << std::endl;
	}

	void setupFromForwardCompact(
		std::string const & rfastaname,
		std::string const & tmpprefix,
		std::string const & rfainame,
		// compact for forward
		std::string const & forwcompact,
		// compact file meta data for forward
		std::string const & forwmetaname,
		// repl fasta for forward
		std::string const & forwreplname,
		int const verbose,
		size_t constrsasamplingrate,
		size_t constrisasamplingrate,
		uint64_t bwtconstrmem,
		unsigned int const numthreads,
		std::ostream & OSI
	)
	{
		fastaname = rfastaname;
		fainame = rfainame;
		metaname = forwmetaname;

		// base name (clip .fasta suffix from index)
		// basefn = libmaus2::util::OutputFileNameTools::clipOff(fastaname,".fasta") + "_rc";
		basefn = tmpprefix + "_rc";
		replname = basefn + ".repl.fasta";
		bpname = basefn + ".repl.fasta.bp";

		histfn = basefn + ".hist";
		isafn = basefn + ".isa";
		preisafn = basefn + ".preisa";
		preisametafn = preisafn + ".meta";

		// generate fasta with replacements
		generateRCFasta(forwreplname, replname);

		compactname = basefn + ".compact";

		loadFAI(fainame,verbose);
		/* static */ generateCompactRCFromCompact(forwcompact, compactname, PFAI.get(), verbose);
		/* static */ generateBP(replname,bpname,verbose);

		//  bwt fn
		bwtfn = basefn + ".bwt";

		generateBwt(compactname,bwtfn,constrsasamplingrate,constrisasamplingrate,bwtconstrmem,numthreads);

		// file name of sampled suffix array (8 bytes per value)
		saname = basefn + ".sa";
		// compacted version of sampled suffix array (minimal number of bits per value in block code)
		compactsaname = saname + ".compact";
		compactsametaname = saname + ".compact.meta";

		n = libmaus2::huffman::RLDecoder::getLength(bwtfn,numthreads);
		nbits = n ? libmaus2::math::numbits(n-1) : 0;
		/* static */ generateCompactSA(saname,compactsaname,compactsametaname,verbose,n,nbits);

		dnarankname = basefn + ".dnarank";
		/* static */ generateDNARank(bwtfn,dnarankname,numthreads,verbose);

		uint64_t s = 0;
		// std::string const combfn = basefn + "_damapper_bwt_forward_index";
		// libmaus2::aio::OutputStreamInstance::unique_ptr_type OSI(new libmaus2::aio::OutputStreamInstance(combfn));
		s += rewriteFAI(OSI,fainame);
		s += rewriteDNARank(OSI,dnarankname);
		s += rewriteDNAIndexMetaDataBigBand(OSI,metaname);
		s += rewriteCompactSampledSuffixArray(OSI,compactsametaname,compactsaname);
		s += rewriteBP(OSI,bpname);
		OSI.flush();
	}

	// main setup for forward index
	void setupFromFasta(
		std::string const & rfastaname,
		std::string const & tmpprefix,
		int const verbose,
		size_t constrsasamplingrate,
		size_t constrisasamplingrate,
		uint64_t bwtconstrmem,
		unsigned int const numthreads,
		std::ostream & OSI
	)
	{
		fastaname = rfastaname;

		fainame = tmpprefix + ".fai";

		// base name (clip .fasta suffix from index)
		basefn = tmpprefix + "_forward"; // libmaus2::util::OutputFileNameTools::clipOff(fastaname,".fasta");

		// FastA with symbol replacements (N->random bases)
		compactname = basefn + ".compact";
		replname = compactname + ".repl.fasta";
		// compact stream meta data name
		metaname = compactname + ".meta";
		//  bwt fn
		bwtfn = basefn + ".bwt";
		// compacted version of replname file
		bpname = compactname + ".repl.fasta.bp";
		// file name of sampled suffix array (8 bytes per value)
		saname = basefn + ".sa";
		// compacted version of sampled suffix array (minimal number of bits per value in block code)
		compactsaname = saname + ".compact";
		compactsametaname = saname + ".compact.meta";
		// dnarank name
		dnarankname = basefn + ".dnarank";

		histfn = basefn + ".hist";
		isafn = basefn + ".isa";
		preisafn = basefn + ".preisa";
		preisametafn = preisafn + ".meta";

		// generate fasta index
		libmaus2::fastx::FastAIndexGenerator::generate(fastaname,fainame,verbose);

		// generate file with base replacements and compact file
		// this generate compactname, replname
		/* static */ generateReplAndCompact(fastaname,replname,compactname);

		// generate BWT + sampled suffix array
		/* static */ generateBwt(compactname,bwtfn,constrsasamplingrate,constrisasamplingrate,bwtconstrmem,numthreads);

		// get length of text
		n = libmaus2::huffman::RLDecoder::getLength(bwtfn,numthreads);
		// number of bits per suffix array element
		nbits = n ? libmaus2::math::numbits(n-1) : 0;

		// generate BP file from replacement fasta
		/* static */ generateBP(replname,bpname,verbose);

		// generate compact suffix array
		/* static */ generateCompactSA(saname,compactsaname,compactsametaname,verbose,n,nbits);

		// generate DNA rank index structure
		generateDNARank(bwtfn,dnarankname,numthreads,verbose);

		uint64_t s = 0;
		// std::string const combfn = basefn + "_damapper_bwt_forward_index";
		// libmaus2::aio::OutputStreamInstance::unique_ptr_type OSI(new libmaus2::aio::OutputStreamInstance(combfn));
		s += rewriteFAI(OSI,fainame);
		s += rewriteDNARank(OSI,dnarankname);
		s += rewriteDNAIndexMetaDataBigBand(OSI,metaname);
		s += rewriteCompactSampledSuffixArray(OSI,compactsametaname,compactsaname);
		s += rewriteBP(OSI,bpname);
		OSI.flush();
	}
};

struct DNAIndex : public DNAIndexBase
{
	std::string membpname;

	libmaus2::fastx::FastAIndex::unique_ptr_type PFAI;
	libmaus2::rank::DNARank::unique_ptr_type PDNA;
	libmaus2::fastx::DNAIndexMetaDataBigBand::unique_ptr_type Pindex;
	libmaus2::rank::DNARankKmerCache::unique_ptr_type PKcache;
	libmaus2::bitio::CompactArray::unique_ptr_type PSSA;
	libmaus2::fastx::CoordinateCache::unique_ptr_type PCC;

	uint64_t n;
	unsigned int nbits;

	uint64_t sasamplingrate;
	uint64_t sasamplingmask;
	unsigned int sashift;

	static std::string getCurrentDir()
	{
		uint64_t size = PATH_MAX;
		while ( true )
		{
			libmaus2::autoarray::AutoArray<char> A(size,false);

			char * c = ::getcwd(A.begin(),A.size());

			if ( ! c )
			{
				int const error = errno;

				switch ( error )
				{
					case ERANGE:
						size *= 2;
						break;
					default:
					{
						libmaus2::exception::LibMausException lme;
						lme.getStream() << "getcwd() failed with error " << strerror(errno) << std::endl;
						lme.finish();
						throw lme;
					}
				}
			}
			else
			{
				return std::string(c);
			}
		}
	}

	static bool startsWithAlphaColon(std::string const & url)
	{
		uint64_t col = url.size();
		for ( uint64_t i = 0; i < url.size() && col == url.size(); ++i )
			if ( url[i] == ':' )
				col = i;

		if ( col == url.size() )
			return false;

		for ( uint64_t i = 0; i < col; ++i )
			if ( !isalpha(static_cast<unsigned char >(url[i])) )
				return false;

		return true;
	}

	libmaus2::bambam::BamHeader::unique_ptr_type getBamHeader(std::map<std::string,std::string> const & m5map, std::string const & refurl, libmaus2::bambam::BamHeader const * inheader, libmaus2::util::ArgParser const & arg)
	{
		std::ostringstream samheaderstr;
		bool havehd = false;

		if ( inheader )
		{
			std::vector<libmaus2::bambam::HeaderLine> VH = libmaus2::bambam::HeaderLine::extractLines(inheader->text);

			for ( uint64_t i = 0; i < VH.size(); ++i )
				if ( VH[i].type != "SQ" )
				{
					samheaderstr << VH[i].line << "\n";
					if ( VH[i].type == "HD" )
						havehd = true;
				}
		}

		if ( ! havehd )
		{
			samheaderstr << "@HD\tVN:1.5\tSO:unknown\n";
		}

		std::ostringstream sqheaderstr;
		std::string storeurl;
		if ( startsWithAlphaColon(refurl) )
			storeurl = refurl;
		else if (refurl.size() && refurl[0] == '/')
			storeurl = std::string("file:") + refurl;
		else
			storeurl = std::string("file:") + getCurrentDir() + "/" + refurl;
		for ( uint64_t i = 0; i < PFAI->size(); ++i )
		{
			libmaus2::fastx::FastAIndexEntry const & entry = (*PFAI)[i];
			sqheaderstr << "@SQ\tSN:" << entry.name << "\tLN:" << entry.length;

			if ( m5map.find(entry.name) != m5map.end() )
				sqheaderstr << "\tM5:" << m5map.find(entry.name)->second;

			sqheaderstr << "\tUR:" << storeurl;

			sqheaderstr << "\n";
		}

		samheaderstr << sqheaderstr.str();

		std::string headertext = samheaderstr.str();

		std::string const upheadtext = ::libmaus2::bambam::ProgramHeaderLineSet::addProgramLine(
			headertext,
			"damapper_bwt", // ID
			"damapper_bwt", // PN
			arg.commandlinecoded, // CL
			::libmaus2::bambam::ProgramHeaderLineSet(headertext).getLastIdInChain(), // PP
			std::string(PACKAGE_VERSION) // VN
		);

		libmaus2::bambam::BamHeader::unique_ptr_type Pbamheader(new libmaus2::bambam::BamHeader(upheadtext));

		Pbamheader->changeSortOrder("unknown");

		return UNIQUE_PTR_MOVE(Pbamheader);
	}

	void loadBP(
		libmaus2::parallel::SimpleThreadPool & STP,
		std::vector<libmaus2::fastx::FastaBPDecoderIdentity::SequenceMeta> & Vseqmeta,
		libmaus2::autoarray::AutoArray<char,libmaus2::autoarray::alloc_type_c> & Aseq,
		uint64_t const numthreads
	)
	{
		libmaus2::aio::InputStreamInstance::unique_ptr_type membpISI(new libmaus2::aio::InputStreamInstance(membpname));
		libmaus2::fastx::FastaBPDecoderIdentity fabpdec(*membpISI);
		fabpdec.decodeSequencesParallel(STP,membpname,numthreads,Aseq,Vseqmeta,true /* map */,4 /* pad symbol */,false /* addrc */);
		assert ( Vseqmeta.size() == PFAI->size() );
	}

	uint64_t getNumPosBits() const
	{
		uint64_t maxl = 0;
		for ( uint64_t i = 0; i < PFAI->size(); ++i )
		{
			uint64_t const l = (*PFAI)[i].length;
			maxl = std::max(l,maxl);
		}
		return maxl ? libmaus2::math::numbits(maxl-1) : 0;
	}

	uint64_t getNumPosBytes() const
	{
		uint64_t const numposbits = getNumPosBits();
		uint64_t const numposbytes = (numposbits+7)/8;
		return numposbytes;
	}

	uint64_t getNumSeqBits() const
	{
		return PFAI->size() ? libmaus2::math::numbits(PFAI->size()-1) : 0;
	}

	uint64_t getNumSeqBytes() const
	{
		uint64_t const numseqbits = getNumSeqBits();
		return (numseqbits+7)/8;
	}

	void fillDazzlerDB(
		HITS_DB & refdb,
		libmaus2::autoarray::AutoArray< HITS_READ > & AREFHITREADS,
		std::vector<libmaus2::fastx::FastaBPDecoderIdentity::SequenceMeta> const & Vseqmeta,
		char * Aseq
	) const
	{
		// set up HITS_DB data structure for reference
		uint64_t refmaxlen = 0;
		uint64_t reftotlen = 0;
		for ( uint64_t i = 0; i < PFAI->size(); ++i )
		{
			uint64_t const l = (*PFAI)[i].length;
			refmaxlen = std::max(refmaxlen,l);
			reftotlen += l;
		}

		AREFHITREADS.ensureSize(PFAI->size() + 2);
		static char refpath[] = "ref.db";

		refdb.ureads = PFAI->size();
		refdb.treads = PFAI->size();
		refdb.cutoff = -1;
		// refdb.all = 0;
		refdb.freq[0] = refdb.freq[1] = refdb.freq[2] = refdb.freq[3] = 0.25;
		refdb.maxlen = refmaxlen;
		refdb.totlen = reftotlen;
		refdb.nreads = PFAI->size();
		refdb.trimmed = 1;
		refdb.part = 0;
		refdb.ufirst = 0;
		refdb.tfirst = 0;
		refdb.path = &refpath[0];
		refdb.loaded = 1;
		refdb.bases = reinterpret_cast<void *>(Aseq);
		refdb.reads = AREFHITREADS.begin() + 1;
		refdb.tracks = NULL;

		for ( uint64_t i = 0; i < Vseqmeta.size(); ++i )
		{
			refdb.reads[i].origin = i;
			refdb.reads[i].rlen = Vseqmeta[i].length;
			refdb.reads[i].fpulse = 0;
			refdb.reads[i].boff = Vseqmeta[i].firstbase;
			assert ( (Aseq + Vseqmeta[i].firstbase)[-1] == 4 );
			refdb.reads[i].coff = -1;
			refdb.reads[i].flags = DB_BEST;
		}

		(reinterpret_cast<int*>(refdb.reads))[-1] = PFAI->size();
		(reinterpret_cast<int*>(refdb.reads))[-2] = PFAI->size();
	}

	void loadFromSingleForward(std::istream & ISI, uint64_t const numthreads, unsigned int const cache_k, int const verbose)
	{
		// FAI
		libmaus2::fastx::FastAIndex::unique_ptr_type TFAI(libmaus2::fastx::FastAIndex::loadSerialised(ISI));
		PFAI = UNIQUE_PTR_MOVE(TFAI);

		// DNARank
		libmaus2::rank::DNARank::unique_ptr_type tPDNA(libmaus2::rank::DNARank::loadFromSerialised(ISI));
		PDNA = UNIQUE_PTR_MOVE(tPDNA);

		// DNAIndexMetaDataBigBand
		libmaus2::fastx::DNAIndexMetaDataBigBand::unique_ptr_type tPindex(libmaus2::fastx::DNAIndexMetaDataBigBand::load(ISI));
		Pindex = UNIQUE_PTR_MOVE(tPindex);

		// CompactArray
		sasamplingrate = libmaus2::util::NumberSerialisation::deserialiseNumber(ISI);
		sasamplingmask = sasamplingrate-1;
		sashift = libmaus2::math::ilog(sasamplingrate);
		assert ( sasamplingrate == 1ull<<sashift );
		libmaus2::bitio::CompactArray::unique_ptr_type tSSA(new libmaus2::bitio::CompactArray(ISI));
		PSSA = UNIQUE_PTR_MOVE(tSSA);

		loadKMerCache(cache_k, numthreads, verbose);
		loadCoordinateCache(verbose);

		uint64_t const bpsize = libmaus2::util::NumberSerialisation::deserialiseNumber(ISI);
		membpname = "mem://damapper_bwt_forw_bp";
		libmaus2::aio::OutputStreamInstance::unique_ptr_type membpOSI(new libmaus2::aio::OutputStreamInstance(membpname));
		libmaus2::util::GetFileSize::copy(ISI,*membpOSI,bpsize);
		membpOSI->flush();
		membpOSI.reset();

		n = PDNA->size();
		nbits = n ? libmaus2::math::numbits(n-1) : 0;
	}

	void loadFromSingleReverseComplement(std::istream & ISI, uint64_t const numthreads, unsigned int const cache_k, int const verbose)
	{
		// FAI
		libmaus2::fastx::FastAIndex::unique_ptr_type TFAI(libmaus2::fastx::FastAIndex::loadSerialised(ISI));
		PFAI = UNIQUE_PTR_MOVE(TFAI);

		// DNARank
		libmaus2::rank::DNARank::unique_ptr_type tPDNA(libmaus2::rank::DNARank::loadFromSerialised(ISI));
		PDNA = UNIQUE_PTR_MOVE(tPDNA);

		// DNAIndexMetaDataBigBand
		libmaus2::fastx::DNAIndexMetaDataBigBand::unique_ptr_type tPindex(libmaus2::fastx::DNAIndexMetaDataBigBand::load(ISI));
		Pindex = UNIQUE_PTR_MOVE(tPindex);

		// CompactArray
		sasamplingrate = libmaus2::util::NumberSerialisation::deserialiseNumber(ISI);
		sasamplingmask = sasamplingrate-1;
		sashift = libmaus2::math::ilog(sasamplingrate);
		assert ( sasamplingrate == 1ull<<sashift );
		libmaus2::bitio::CompactArray::unique_ptr_type tSSA(new libmaus2::bitio::CompactArray(ISI));
		PSSA = UNIQUE_PTR_MOVE(tSSA);

		loadKMerCache(cache_k, numthreads, verbose);
		loadCoordinateCache(verbose);

		uint64_t const bpsize = libmaus2::util::NumberSerialisation::deserialiseNumber(ISI);
		membpname = "mem://damapper_bwt_reco_bp";
		libmaus2::aio::OutputStreamInstance::unique_ptr_type membpOSI(new libmaus2::aio::OutputStreamInstance(membpname));
		libmaus2::util::GetFileSize::copy(ISI,*membpOSI,bpsize);
		membpOSI->flush();
		membpOSI.reset();

		n = PDNA->size();
		nbits = n ? libmaus2::math::numbits(n-1) : 0;
	}


	DNAIndex()
	{

	}

	uint64_t lookupSA(uint64_t const jr) const
	{
		if ( sasamplingrate == 1 )
		{
			return (*(PSSA))[jr];
		}
		else
		{
			std::pair<uint64_t,uint64_t> const SP = PDNA->simpleLFUntilMask(jr, sasamplingmask);
			return ((*(PSSA))[SP.first >> sashift] + SP.second)%PDNA->size();
		}
	}

	void loadKMerCache(unsigned int const cache_k, uint64_t const numthreads, int const verbose)
	{
		if ( verbose > 1 )
			std::cerr << "[V] computing " << cache_k << "-mer cache...";
		libmaus2::rank::DNARankKmerCache::unique_ptr_type tKcache(new libmaus2::rank::DNARankKmerCache(*PDNA,cache_k,numthreads));
		PKcache = UNIQUE_PTR_MOVE(tKcache);
		if ( verbose > 1 )
			std::cerr << "done." << std::endl;
	}

	void loadCoordinateCache(int const verbose)
	{
		if ( verbose > 1 )
			std::cerr << "[V] constructing coordinate cache...";
		libmaus2::fastx::CoordinateCache::unique_ptr_type tCC(new libmaus2::fastx::CoordinateCache(*PDNA,*Pindex,16));
		PCC = UNIQUE_PTR_MOVE(tCC);
		if ( verbose > 1 )
			std::cerr << "done." << std::endl;
	}
};

struct ReadInterval
{
	uint64_t from;
	uint64_t to;
	uint64_t id;

	ReadInterval() {}
	ReadInterval(uint64_t const rfrom, uint64_t const rto, uint64_t const rid) : from(rfrom), to(rto), id(rid) {}

	uint64_t getFrom() const
	{
		return from;
	}

	uint64_t getTo() const
	{
		return to;
	}
};

struct RefKMer
{
	uint64_t code;
	uint64_t low;
	uint64_t high;
};

struct RefKMerCodeComparator
{
	bool operator()(RefKMer const & A, RefKMer const & B) const
	{
		return A.code < B.code;
	}
};

struct QueryKMer
{
	uint64_t code;
	uint32_t rpos;
	uint32_t readid;

	bool operator<(QueryKMer const & O) const
	{
		if ( code != O.code )
			return code < O.code;
		else if ( readid != O.readid )
			return readid < O.readid;
		else
			return rpos < O.rpos;
	}
};

struct PackageSearchKmers
{
	uint64_t t; /* thread id */
	uint64_t readsperpack;
	uint64_t numreads;
	unsigned int k;
	libmaus2::autoarray::AutoArray< libmaus2::autoarray::AutoArray < RefKMer >::unique_ptr_type > * RKM;
	libmaus2::autoarray::AutoArray< libmaus2::autoarray::AutoArray < QueryKMer >::unique_ptr_type > * QKM;
	libmaus2::autoarray::AutoArray<char> const * Areaddata;
	libmaus2::autoarray::AutoArray< libmaus2::fastx::LineBufferFastAReader::ReadMeta > const * O;
	libmaus2::autoarray::AutoArray< uint64_t > * NUMRK;
	libmaus2::autoarray::AutoArray< uint64_t > * NUMQK;
	DNAIndex * index;
	libmaus2::parallel::PosixSemaphore *finsem;

	PackageSearchKmers() {}

	PackageSearchKmers(
		uint64_t const rt, /* thread id */
		uint64_t const rreadsperpack,
		uint64_t const rnumreads,
		unsigned int const rk,
		libmaus2::autoarray::AutoArray< libmaus2::autoarray::AutoArray < RefKMer >::unique_ptr_type > * rRKM,
		libmaus2::autoarray::AutoArray< libmaus2::autoarray::AutoArray < QueryKMer >::unique_ptr_type > * rQKM,
		libmaus2::autoarray::AutoArray<char> const * rAreaddata,
		libmaus2::autoarray::AutoArray< libmaus2::fastx::LineBufferFastAReader::ReadMeta > const * rO,
		libmaus2::autoarray::AutoArray< uint64_t > * rNUMRK,
		libmaus2::autoarray::AutoArray< uint64_t > * rNUMQK,
		DNAIndex * rindex,
		libmaus2::parallel::PosixSemaphore *rfinsem
	)
	: t(rt), readsperpack(rreadsperpack), numreads(rnumreads), k(rk), RKM(rRKM), QKM(rQKM), Areaddata(rAreaddata), O(rO), NUMRK(rNUMRK), NUMQK(rNUMQK), index(rindex), finsem(rfinsem)
	{

	}
};

struct PackageCompressBgzf
{
	libmaus2::parallel::LockedGrowingFreeList<libmaus2::lz::BgzfDeflateZStreamBase,BgzfDeflateZStreamBaseAllocator,BgzfDeflateZStreamBaseTypeInfo> * zstreambasefreelist;
	std::pair<uint8_t *,uint8_t *> P;
	libmaus2::lz::BgzfDeflateOutputBufferBase::shared_ptr_type B;
	libmaus2::lz::BgzfDeflateZStreamBaseFlushInfo * info;
	libmaus2::parallel::PosixSemaphore *finsem;

	PackageCompressBgzf() {}
	PackageCompressBgzf(
		libmaus2::parallel::LockedGrowingFreeList<libmaus2::lz::BgzfDeflateZStreamBase,BgzfDeflateZStreamBaseAllocator,BgzfDeflateZStreamBaseTypeInfo> * rzstreambasefreelist,
		std::pair<uint8_t *,uint8_t *> rP,
		libmaus2::lz::BgzfDeflateOutputBufferBase::shared_ptr_type rB,
		libmaus2::lz::BgzfDeflateZStreamBaseFlushInfo * rinfo,
		libmaus2::parallel::PosixSemaphore *rfinsem
	) : zstreambasefreelist(rzstreambasefreelist), P(rP), B(rB), info(rinfo), finsem(rfinsem)
	{}

	void dispatch()
	{
		libmaus2::lz::BgzfDeflateZStreamBase::shared_ptr_type zstreambase = zstreambasefreelist->get();
		*info = zstreambase->flush(P.first,P.second,*B);
		zstreambasefreelist->put(zstreambase);
		finsem->post();
	}
};

static void threadSearchKmers(PackageSearchKmers package)
{
	uint64_t const low = package.t * package.readsperpack;
	uint64_t const high = std::min(low+package.readsperpack,package.numreads);
	libmaus2::autoarray::AutoArray < RefKMer > & LRKM = *(*(package.RKM))[package.t];
	libmaus2::autoarray::AutoArray < QueryKMer > & LQKM = *(*(package.QKM))[package.t];

	uint64_t rkmo = 0;
	uint64_t qkmo = 0;

	uint64_t const kmask = libmaus2::math::lowbits(2*package.k);

	for ( uint64_t zz = low; zz < high; ++zz )
	{
		//char const * const cmapped = Areaddata.begin() + O[zz].dataoff;
		char const * const cmapped = package.Areaddata->begin() + (*(package.O))[zz].dataoff;
		uint64_t const m = (*(package.O))[zz].len;

		uint64_t const numk = (m >= package.k) ? m-package.k+1 : 0;

		if ( numk && package.k )
		{
			uint64_t v = 0;
			uint64_t e = 0;
			char const * c = cmapped;

			for ( uint64_t i = 0; i < package.k-1; ++i )
			{
				v <<= 2;
				uint64_t const s = *(c++);
				v |= s & 3;
				e += (s>>2)&1;
			}
			char const * b = cmapped;

			for ( uint64_t ik = 0; ik < numk; ++ik )
			{
				v <<= 2;
				v &= kmask;
				uint64_t const s = *(c++);
				v |= s & 3;
				e += (s >> 2)&1;

				if ( ! e )
				{
					std::pair<uint64_t,uint64_t> const P = package.index->PKcache->search(c-package.k,package.k);

					if ( P.second != P.first )
					{
						RefKMer K;
						K.code = v;
						K.low = P.first;
						K.high = P.second;
						LRKM.push(rkmo,K);

						assert ( (K.code & kmask) == K.code );

						QueryKMer Q;
						Q.code = v;
						Q.rpos = ik;
						Q.readid = zz;
						LQKM.push(qkmo,Q);
					}
				}

				uint64_t sa = *(b++);
				e -= (sa >> 2)&1;
			}
		}
	}

	// sort reference kmers
	std::sort(LRKM.begin(),LRKM.begin() + rkmo,RefKMerCodeComparator());
	uint64_t l = 0;
	uint64_t o = 0;
	// uniquify
	while ( l < rkmo )
	{
		uint64_t h = l+1;
		while ( h < rkmo && LRKM[h].code == LRKM[l].code )
			++h;
		LRKM[o++] = LRKM[l];
		l = h;
	}
	rkmo = o;

	(*(package.NUMRK))[package.t] = rkmo;
	(*(package.NUMQK))[package.t] = qkmo;

	package.finsem->post();
}

struct PackageConcatSA
{
	uint64_t t;
	libmaus2::autoarray::AutoArray < RefKMer > * GRKM;
	libmaus2::autoarray::AutoArray< uint64_t > * NUMRK;
	libmaus2::autoarray::AutoArray< libmaus2::autoarray::AutoArray < RefKMer >::unique_ptr_type > * RKM;
	libmaus2::parallel::PosixSemaphore *finsem;

	PackageConcatSA() {}

	PackageConcatSA(uint64_t const rt,
		libmaus2::autoarray::AutoArray < RefKMer > * rGRKM,
		libmaus2::autoarray::AutoArray< uint64_t > * rNUMRK,
		libmaus2::autoarray::AutoArray< libmaus2::autoarray::AutoArray < RefKMer >::unique_ptr_type > * rRKM,
		libmaus2::parallel::PosixSemaphore *rfinsem)
	: t(rt), GRKM(rGRKM), NUMRK(rNUMRK), RKM(rRKM), finsem(rfinsem) {}

};

static void threadConcatSA(PackageConcatSA package)
{
	RefKMer * O = (package.GRKM)->begin() + (*(package.NUMRK))[(package.t)];
	RefKMer * Oe = (package.GRKM)->begin() + (*(package.NUMRK))[(package.t)+1];
	RefKMer * I = (*(package.RKM))[(package.t)]->begin();

	while ( O != Oe )
		*(O++) = *(I++);

	(*(package.RKM))[(package.t)]->resize(0);

	package.finsem->post();
}

struct PackageFillHist
{
	uint64_t t;
	libmaus2::autoarray::AutoArray< libmaus2::autoarray::AutoArray<uint64_t>::unique_ptr_type> * Ahist;
	libmaus2::parallel::PosixSemaphore *finsem;

	PackageFillHist() {}

	PackageFillHist(
		uint64_t const rt,
		libmaus2::autoarray::AutoArray< libmaus2::autoarray::AutoArray<uint64_t>::unique_ptr_type> * rAhist,
		libmaus2::parallel::PosixSemaphore *rfinsem
	) : t(rt), Ahist(rAhist), finsem(rfinsem) {}
};

static void threadFillHist(PackageFillHist package)
{
	std::fill((*(package.Ahist))[(package.t)]->begin(),(*(package.Ahist))[(package.t)]->end(),0ull);
	package.finsem->post();
}

struct PackageUnifyPackageSizes
{
	uint64_t t;
	libmaus2::autoarray::AutoArray< libmaus2::autoarray::AutoArray<uint64_t>::unique_ptr_type> * Ahist;
	libmaus2::autoarray::AutoArray < RefKMer > * GRKM;
	std::vector<uint64_t> * ubounds;
	std::vector<uint64_t> * usize;
	uint64_t histthres;
	libmaus2::parallel::PosixSemaphore *finsem;

	PackageUnifyPackageSizes() {}

	PackageUnifyPackageSizes(
		uint64_t const rt,
		libmaus2::autoarray::AutoArray< libmaus2::autoarray::AutoArray<uint64_t>::unique_ptr_type> * rAhist,
		libmaus2::autoarray::AutoArray < RefKMer > * rGRKM,
		std::vector<uint64_t> * rubounds,
		std::vector<uint64_t> * rusize,
		uint64_t const rhistthres,
		libmaus2::parallel::PosixSemaphore *rfinsem
	) : t(rt), Ahist(rAhist), GRKM(rGRKM), ubounds(rubounds), usize(rusize), histthres(rhistthres), finsem(rfinsem) {}
};

static void threadUnifyPackageSizes(PackageUnifyPackageSizes package)
{
	libmaus2::autoarray::AutoArray<uint64_t> & thist = *((*(package.Ahist))[(package.t)-1].get());

	// lower bound
	uint64_t ulow = (*(package.ubounds))[(package.t)-1];
	// upper bound
	uint64_t const utop = (*(package.ubounds))[(package.t)];

	// output pointer
	uint64_t uout = ulow;
	while ( ulow < utop )
	{
		// go to end of code
		uint64_t uhigh = ulow+1;
		while ( uhigh < utop && (*(package.GRKM))[uhigh].code == (*(package.GRKM))[ulow].code )
			++uhigh;

		// single output
		(*(package.GRKM))[uout++] = (*(package.GRKM))[ulow];

		uint64_t const s = (*(package.GRKM))[ulow].high-(*(package.GRKM))[ulow].low;
		if ( s < (package.histthres) )
			thist [ s ] += 1;
		else
			thist [ (package.histthres) ] += s;

		ulow = uhigh;
	}

	// size of package
	(*(package.usize))[(package.t)-1] = uout - (*(package.ubounds))[(package.t)-1];

	package.finsem->post();
}

struct PackageMergeHistograms
{
	uint64_t t;
	uint64_t tperthread;
	uint64_t histmult;
	uint64_t numthreads;
	libmaus2::autoarray::AutoArray< libmaus2::autoarray::AutoArray<uint64_t>::unique_ptr_type> * Ahist;
	libmaus2::parallel::PosixSemaphore *finsem;

	PackageMergeHistograms() {}

	PackageMergeHistograms(
		uint64_t const rt,
		uint64_t const rtperthread,
		uint64_t const rhistmult,
		uint64_t const rnumthreads,
		libmaus2::autoarray::AutoArray< libmaus2::autoarray::AutoArray<uint64_t>::unique_ptr_type> * rAhist,
		libmaus2::parallel::PosixSemaphore *rfinsem
	) : t(rt), tperthread(rtperthread), histmult(rhistmult), numthreads(rnumthreads), Ahist(rAhist), finsem(rfinsem)
	{}
};

static void threadMergeHistograms(PackageMergeHistograms package)
{
	uint64_t const tlow = package.t * package.tperthread;
	uint64_t const thigh = std::min ( tlow + package.tperthread, package.histmult );

	uint64_t * const T = (*(package.Ahist))[0]->begin();

	for ( uint64_t j = 1; j < package.numthreads; ++j )
	{
		uint64_t * const S = (*(package.Ahist))[j]->begin();

		for ( uint64_t i = tlow; i < thigh; ++i )
			T[i] += S[i];
	}

	package.finsem->post();
}

struct PackageConcatenateUnique
{
	uint64_t t;
	libmaus2::autoarray::AutoArray < RefKMer > * GRKM;
	libmaus2::autoarray::AutoArray < RefKMer > * TGRKM;
	std::vector<uint64_t> * ubounds;
	std::vector<uint64_t> * usize;
	libmaus2::parallel::PosixSpinLock * gclock;
	uint64_t volatile * gc;
	libmaus2::parallel::PosixSemaphore *finsem;

	PackageConcatenateUnique() {}
	PackageConcatenateUnique(
		uint64_t const rt,
		libmaus2::autoarray::AutoArray < RefKMer > * rGRKM,
		libmaus2::autoarray::AutoArray < RefKMer > * rTGRKM,
		std::vector<uint64_t> * rubounds,
		std::vector<uint64_t> * rusize,
		libmaus2::parallel::PosixSpinLock * rgclock,
		uint64_t volatile * rgc,
		libmaus2::parallel::PosixSemaphore *rfinsem
	) : t(rt), GRKM(rGRKM), TGRKM(rTGRKM), ubounds(rubounds), usize(rusize), gclock(rgclock), gc(rgc), finsem(rfinsem) {}
};

static void threadConcatenateUnique(
	PackageConcatenateUnique package
	#if 0
	uint64_t const t,
	libmaus2::autoarray::AutoArray < RefKMer > * GRKM,
	libmaus2::autoarray::AutoArray < RefKMer > * TGRKM,
	std::vector<uint64_t> * ubounds,
	std::vector<uint64_t> * usize,
	libmaus2::parallel::PosixSpinLock * gclock,
	uint64_t volatile * gc,
	libmaus2::parallel::PosixSemaphore *finsem
	#endif
)
{
	RefKMer * I = (package.GRKM)->begin() + (*(package.ubounds))[(package.t)-1];
	RefKMer * O = (package.TGRKM)->begin() + (*(package.usize))[(package.t)-1];
	RefKMer * Oe = (package.TGRKM)->begin() + (*(package.usize))[(package.t)];

	uint64_t c = 0;
	while ( O != Oe )
	{
		c += I->high-I->low;
		*(O++) = *(I++);
	}
	package.gclock->lock();
	(*(package.gc)) += c;
	package.gclock->unlock();

	#if 0
	uint64_t cc = 0;
	for ( RefKMer * C = TGRKM.begin() + (package.usize)[t-1]; C != TGRKM.begin() + (package.usize)[t]; ++C )
		cc += C->high-C->low;
	assert ( cc == c );
	#endif

	package.finsem->post();
}

struct PackageSumIntervals
{
	uint64_t t;
	uint64_t supacksize;
	std::vector<uint64_t> * usize;
	libmaus2::autoarray::AutoArray < RefKMer > * GRKM;
	libmaus2::autoarray::AutoArray<uint64_t,libmaus2::autoarray::alloc_type_c> * PP;
	libmaus2::parallel::PosixSemaphore *finsem;

	PackageSumIntervals() {}
	PackageSumIntervals(
		uint64_t const rt,
		uint64_t const rsupacksize,
		std::vector<uint64_t> * rusize,
		libmaus2::autoarray::AutoArray < RefKMer > * rGRKM,
		libmaus2::autoarray::AutoArray<uint64_t,libmaus2::autoarray::alloc_type_c> * rPP,
		libmaus2::parallel::PosixSemaphore *rfinsem
	) : t(rt), supacksize(rsupacksize), usize(rusize), GRKM(rGRKM), PP(rPP), finsem(rfinsem)
	{

	}
};

static void threadSumIntervals(
	PackageSumIntervals package
	#if 0
	uint64_t const t,
	uint64_t const supacksize,
	std::vector<uint64_t> * usize,
	libmaus2::autoarray::AutoArray < RefKMer > * GRKM,
	libmaus2::autoarray::AutoArray<uint64_t,libmaus2::autoarray::alloc_type_c> * PP,
	libmaus2::parallel::PosixSemaphore *finsem
	#endif
)
{
	uint64_t const ulow = package.t * package.supacksize;
	uint64_t const uhigh = std::min(ulow+package.supacksize,package.usize->back());
	assert ( uhigh != ulow );

	RefKMer * I  = package.GRKM->begin() + ulow;
	RefKMer * Ie = package.GRKM->begin() + uhigh;

	uint64_t s = 0;
	for ( ; I != Ie ; ++I )
		s += I->high-I->low;
	(*(package.PP))[uhigh] = s;

	package.finsem->post();
}

struct PackagePrefixSums
{
	uint64_t t;
	uint64_t supacksize;
	std::vector<uint64_t> * usize;
	libmaus2::autoarray::AutoArray < RefKMer > * GRKM;
	libmaus2::autoarray::AutoArray<uint64_t,libmaus2::autoarray::alloc_type_c> * PP;
	libmaus2::parallel::PosixSemaphore *finsem;

	PackagePrefixSums() {}

	PackagePrefixSums(
		uint64_t const rt,
		uint64_t const rsupacksize,
		std::vector<uint64_t> * rusize,
		libmaus2::autoarray::AutoArray < RefKMer > * rGRKM,
		libmaus2::autoarray::AutoArray<uint64_t,libmaus2::autoarray::alloc_type_c> * rPP,
		libmaus2::parallel::PosixSemaphore *rfinsem
	) : t(rt), supacksize(rsupacksize), usize(rusize), GRKM(rGRKM), PP(rPP), finsem(rfinsem)
	{

	}
};

static void threadPrefixSums(PackagePrefixSums package)
{
	uint64_t const ulow = package.t * package.supacksize;
	uint64_t const uhigh = std::min(ulow+package.supacksize,package.usize->back());
	assert ( uhigh != ulow );

	RefKMer * I  = package.GRKM->begin() + ulow;
	RefKMer * Ie = package.GRKM->begin() + uhigh;

	uint64_t * p = package.PP->begin() + ulow;

	uint64_t s = *p;
	for ( ; I != Ie ; ++I )
	{
		*(p++) = s;
		s += (I->high-I->low);
	}

	package.finsem->post();
}

struct PackageLookupSA
{
	uint64_t t;
	unsigned int k;
	libmaus2::autoarray::AutoArray < RefKMer > * GRKM;
	libmaus2::autoarray::AutoArray < QueryKMer, libmaus2::autoarray::alloc_type_c > * GRQKM;
	std::vector<uint64_t> * uranges;
	std::vector<uint64_t> * urefk;
	std::vector<uint64_t> * NUMREFKVALID;
	DNAIndex * index;
	libmaus2::parallel::PosixSemaphore *finsem;

	PackageLookupSA() {}

	PackageLookupSA(
		uint64_t const rt,
		unsigned int const rk,
		libmaus2::autoarray::AutoArray < RefKMer > * rGRKM,
		libmaus2::autoarray::AutoArray < QueryKMer, libmaus2::autoarray::alloc_type_c > * rGRQKM,
		std::vector<uint64_t> * ruranges,
		std::vector<uint64_t> * rurefk,
		std::vector<uint64_t> * rNUMREFKVALID,
		DNAIndex * rindex,
		libmaus2::parallel::PosixSemaphore *rfinsem
	) : t(rt), k(rk), GRKM(rGRKM), GRQKM(rGRQKM), uranges(ruranges), urefk(rurefk), NUMREFKVALID(rNUMREFKVALID), index(rindex), finsem(rfinsem)
	{}

};

static void threadLookupSA(PackageLookupSA package)
{
	RefKMer * I  = package.GRKM->begin() + (*(package.uranges))[package.t-1];
	RefKMer * Ie = package.GRKM->begin() + (*(package.uranges))[package.t];
	QueryKMer * O = package.GRQKM->begin() + (*(package.urefk))[package.t-1];
	QueryKMer Q;
	uint64_t numvalid = 0;
	uint64_t c = 0;

	for ( ; I != Ie; ++I )
	{
		RefKMer const & R = *I;
		c += R.high-R.low;
		Q.code = R.code;

		for ( uint64_t jr = R.low; jr < R.high; ++jr )
		{
			uint64_t const p = package.index->lookupSA(jr);
			std::pair<uint64_t,uint64_t> coord = (*(package.index->PCC))[p]; // Pindex->mapCoordinates((*SSA)[jr]);
			if ( package.index->Pindex->valid(coord,package.k) )
			{
				Q.rpos = coord.second;
				Q.readid = coord.first;
				*(O++) = Q;
				numvalid += 1;
			}
			#if 0
			else
			{
				Q.code = std::numeric_limits<uint64_t>::max();
				Q.rpos = 0;
				Q.readid = 0;
				*(O++) = Q;
			}
			#endif
		}
	}
	assert ( c == (*(package.urefk))[package.t]-(*(package.urefk))[package.t-1] );
	// assert ( O == GRQKM.begin() + NUMREFK[t] );
	(*(package.NUMREFKVALID))[package.t-1] = numvalid;

	package.finsem->post();
}

struct PackageCompactPos
{
	uint64_t t;
	libmaus2::autoarray::AutoArray < QueryKMer, libmaus2::autoarray::alloc_type_c > * GRQKM;
	libmaus2::autoarray::AutoArray < QueryKMer, libmaus2::autoarray::alloc_type_c > * TGRQKM;
	std::vector<uint64_t> * NUMREFKVALID;
	std::vector<uint64_t> * urefk;
	libmaus2::parallel::PosixSemaphore *finsem;

	PackageCompactPos() {}
	PackageCompactPos(
		uint64_t const rt,
		libmaus2::autoarray::AutoArray < QueryKMer, libmaus2::autoarray::alloc_type_c > * rGRQKM,
		libmaus2::autoarray::AutoArray < QueryKMer, libmaus2::autoarray::alloc_type_c > * rTGRQKM,
		std::vector<uint64_t> * rNUMREFKVALID,
		std::vector<uint64_t> * rurefk,
		libmaus2::parallel::PosixSemaphore *rfinsem
	) : t(rt), GRQKM(rGRQKM), TGRQKM(rTGRQKM), NUMREFKVALID(rNUMREFKVALID), urefk(rurefk), finsem(rfinsem)
	{}
};

static void threadCompactPos(PackageCompactPos package)
{
	QueryKMer * I = package.GRQKM->begin() + (*(package.urefk))[package.t-1];
	QueryKMer * O = package.TGRQKM->begin() + (*(package.NUMREFKVALID))[package.t-1];
	QueryKMer * Oe = package.TGRQKM->begin() + (*(package.NUMREFKVALID))[package.t];

	while ( O != Oe )
		*(O++) = *(I++);

	package.finsem->post();
}

struct PackageConcatQueryKmers
{
	uint64_t t;
	libmaus2::autoarray::AutoArray< libmaus2::autoarray::AutoArray < QueryKMer >::unique_ptr_type > * QKM;
	libmaus2::autoarray::AutoArray < QueryKMer, libmaus2::autoarray::alloc_type_c > * GQKM;
	libmaus2::autoarray::AutoArray< uint64_t > * NUMQK;
	libmaus2::parallel::PosixSemaphore *finsem;

	PackageConcatQueryKmers() {}
	PackageConcatQueryKmers(
		uint64_t const rt,
		libmaus2::autoarray::AutoArray< libmaus2::autoarray::AutoArray < QueryKMer >::unique_ptr_type > * rQKM,
		libmaus2::autoarray::AutoArray < QueryKMer, libmaus2::autoarray::alloc_type_c > * rGQKM,
		libmaus2::autoarray::AutoArray< uint64_t > * rNUMQK,
		libmaus2::parallel::PosixSemaphore *rfinsem
	) : t(rt), QKM(rQKM), GQKM(rGQKM), NUMQK(rNUMQK), finsem(rfinsem)
	{}
};

static void threadConcatQueryKmers(PackageConcatQueryKmers package)
{
	libmaus2::autoarray::AutoArray < QueryKMer > & LQKM = *(*(package.QKM))[package.t];

	QueryKMer const * I = LQKM.begin();
	QueryKMer * O = package.GQKM->begin() + (*(package.NUMQK))[package.t];
	QueryKMer * Oe = package.GQKM->begin() + (*(package.NUMQK))[package.t+1];

	while ( O != Oe )
		*(O++) = *(I++);

	(*(package.QKM))[package.t]->resize(0);

	package.finsem->post();
}

struct LASToBamContext
{
	typedef LASToBamContext this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;

	libmaus2::dazzler::align::LASToBamConverterBase converter;
	libmaus2::bambam::parallel::FragmentAlignmentBufferFragment fragment;
	libmaus2::autoarray::AutoArray<char> ARC;

	LASToBamContext(
		int64_t const rtspace,
		bool const rcalmdnm,
		libmaus2::dazzler::align::LASToBamConverterBase::supplementary_seq_strategy_t rseqstrat,
		std::string const & rrgid,
		libmaus2::dazzler::align::RefMapEntryVector const & refmap
	) : converter(rtspace,rcalmdnm,rseqstrat,rrgid,refmap), fragment(), ARC()
	{

	}
};

struct PackageEncodeBam
{
	uint64_t zr;
	std::vector<uint64_t> const * ointv;
	libmaus2::autoarray::AutoArray<LASToBamContext::unique_ptr_type> * Acontexts;
	libmaus2::dazzler::align::OverlapData const * OVLdata;
	int64_t prevmark;
	libmaus2::autoarray::AutoArray<char> const * Areadnames;
	libmaus2::bambam::BamHeader::unique_ptr_type * Pbamheader;
	libmaus2::autoarray::AutoArray< libmaus2::fastx::LineBufferFastAReader::ReadMeta > const * O;
	HITS_DB * refdb;
	HITS_DB * readsdb;

	libmaus2::autoarray::AutoArray<uint8_t> * Abam;
	libmaus2::autoarray::AutoArray< uint64_t > * Obam;
	uint64_t oobam;

	libmaus2::parallel::PosixSemaphore *finsem;

	PackageEncodeBam() {}

	PackageEncodeBam(
		uint64_t const rzr,
		std::vector<uint64_t> const * rointv,
		libmaus2::autoarray::AutoArray<LASToBamContext::unique_ptr_type> * rAcontexts,
		libmaus2::dazzler::align::OverlapData const * rOVLdata,
		int64_t const rprevmark,
		libmaus2::autoarray::AutoArray<char> const * rAreadnames,
		libmaus2::bambam::BamHeader::unique_ptr_type * rPbamheader,
		libmaus2::autoarray::AutoArray< libmaus2::fastx::LineBufferFastAReader::ReadMeta > const * rO,
		HITS_DB * rrefdb,
		HITS_DB * rreadsdb,
		libmaus2::autoarray::AutoArray<uint8_t> * rAbam,
		libmaus2::autoarray::AutoArray< uint64_t > * rObam,
		uint64_t roobam,
		libmaus2::parallel::PosixSemaphore *rfinsem
	) : zr(rzr), ointv(rointv), Acontexts(rAcontexts), OVLdata(rOVLdata), prevmark(rprevmark), Areadnames(rAreadnames),
	    Pbamheader(rPbamheader), O(rO), refdb(rrefdb), readsdb(rreadsdb), Abam(rAbam), Obam(rObam), oobam(roobam), finsem(rfinsem)
	{

	}
};

static void threadEncodeBam(PackageEncodeBam package)
{
	uint64_t tid = package.zr-1;
	uint64_t const rlow = (*(package.ointv))[package.zr-1];
	uint64_t const rhigh = (*(package.ointv))[package.zr];
	LASToBamContext & context = *((*(package.Acontexts))[tid]);
	context.fragment.reset();

	// last read id for previous block
	int64_t const prevb = rlow ? static_cast<int64_t>(libmaus2::dazzler::align::OverlapData::getBRead((package.OVLdata)->getData(rlow-1).first)) : package.prevmark;
	// next expected read
	int64_t nextexpt = prevb+1;

	//std::cerr << "tid " << tid << " " << nextexpt << std::endl;
	libmaus2::bambam::BamAuxFilterVector auxfilter;
	auxfilter.set('A','S');
	auxfilter.set('M','D');
	auxfilter.set('N','M');
	auxfilter.set('c','i');
	auxfilter.set('c','n');
	auxfilter.set('c','j');
	auxfilter.set('c','l');
	libmaus2::autoarray::AutoArray < libmaus2::bambam::BamAlignmentDecoderBase::AuxInfo > auxinfo;
	libmaus2::autoarray::AutoArray < libmaus2::dazzler::align::LASToBamConverterBase::AuxTagAddRequest const * > Aauxadd;
	libmaus2::autoarray::AutoArray < libmaus2::dazzler::align::LASToBamConverterBase::AuxTagCopyAddRequest > Aauxcopy;

	// process in batches of equal b read id
	uint64_t ilow = rlow;
	while ( ilow < rhigh )
	{
		uint64_t ihigh = ilow+1;
		int64_t const refbread = libmaus2::dazzler::align::OverlapData::getBRead((package.OVLdata)->getData(ilow).first);
		while ( ihigh < rhigh && libmaus2::dazzler::align::OverlapData::getBRead((package.OVLdata)->getData(ihigh).first) == refbread )
			++ihigh;

		while ( nextexpt < refbread )
		{
			char const * rname = (package.Areadnames)->begin() + (*(package.O))[nextexpt].nameoff;
			char const * rdata = reinterpret_cast<char const *>((package.readsdb)->bases) + (package.readsdb)->reads[nextexpt].boff;
			uint64_t const readlen = (package.readsdb)->reads[nextexpt].rlen;
			assert ( rdata[-1] == 4 );
			assert ( rdata[readlen] == 4 );

			uint8_t * pbam = 0;
			uint64_t bamlen = 0;
			uint64_t oauxadd = 0;
			uint64_t oauxcopy = 0;

			if ( nextexpt + 1 < static_cast<int64_t>(package.oobam) )
			{
				pbam = package.Abam->begin() + package.Obam->at(nextexpt);
				bamlen = package.Obam->at(nextexpt+1)-package.Obam->at(nextexpt);

				uint64_t const naux = libmaus2::bambam::BamAlignmentDecoderBase::enumerateAuxTagsFilterOut(pbam,bamlen,auxinfo,auxfilter);

				for ( uint64_t i = 0; i < naux; ++i )
				{
					Aauxcopy.push(
						oauxcopy,
						libmaus2::dazzler::align::LASToBamConverterBase::AuxTagCopyAddRequest(pbam + auxinfo[i].o,auxinfo[i].l)
					);
					// std::cerr << "[" << i << "] " << auxinfo[i].tag[0] << auxinfo[i].tag[1] << std::endl;
				}
			}

			for ( uint64_t i = 0; i < oauxcopy; ++i )
				Aauxadd.push(oauxadd,&Aauxcopy[i]);

			libmaus2::dazzler::align::LASToBamConverterBase::AuxTagAddRequest const ** aux_a = &Aauxadd[0];
			libmaus2::dazzler::align::LASToBamConverterBase::AuxTagAddRequest const ** aux_e = &Aauxadd[oauxadd];

			context.converter.convertUnmapped(rdata,readlen,rname,context.fragment,**(package.Pbamheader),aux_a,aux_e);

			// std::cerr << "[D] read " << nextexpt << " " << rname << " was not aligned" << " ref " << refbread << std::endl;

			++nextexpt;
		}

		assert ( nextexpt == refbread );

		uint64_t numchains = 0;
		for ( uint64_t z = ilow; z < ihigh; ++z )
		{
			std::pair<uint8_t const *, uint8_t const *> const P = (package.OVLdata)->getData(z);
			bool const isStart = libmaus2::dazzler::align::OverlapData::getStartFlag(P.first);
			if ( isStart )
				numchains++;
		}

		uint8_t * pbam = 0;
		uint64_t bamlen = 0;

		if ( refbread + 1 < static_cast<int64_t>(package.oobam) )
		{
			pbam = package.Abam->begin() + package.Obam->at(refbread);
			bamlen = package.Obam->at(refbread+1)-package.Obam->at(refbread);
		}

		uint64_t il = ilow;
		uint64_t chainid = 0;
		while ( il < ihigh )
		{
			uint64_t ih = il+1;
			while (
				ih < ihigh
				&&
				(!libmaus2::dazzler::align::OverlapData::getStartFlag((package.OVLdata)->getData(ih).first))
			)
				++ih;

			uint64_t const lchainid = chainid++;
			bool const secondary = (lchainid > 0);

			for ( uint64_t z = il; z < ih; ++z )
			{
				std::pair<uint8_t const *, uint8_t const *> const P = (package.OVLdata)->getData(z);

				int64_t const aread = libmaus2::dazzler::align::OverlapData::getARead(P.first);
				int64_t const bread = libmaus2::dazzler::align::OverlapData::getBRead(P.first);
				int64_t const flags = libmaus2::dazzler::align::OverlapData::getFlags(P.first);
				bool const inverse = libmaus2::dazzler::align::OverlapData::getInverseFlag(P.first);
				char const * readname = (package.Areadnames)->begin() + (*(package.O))[bread].nameoff;

				uint64_t oauxcopy = 0;
				if ( pbam )
				{
					uint64_t const naux = libmaus2::bambam::BamAlignmentDecoderBase::enumerateAuxTagsFilterOut(pbam,bamlen,auxinfo,auxfilter);

					for ( uint64_t i = 0; i < naux; ++i )
					{
						Aauxcopy.push(
							oauxcopy,
							libmaus2::dazzler::align::LASToBamConverterBase::AuxTagCopyAddRequest(pbam + auxinfo[i].o,auxinfo[i].l)
						);
						// std::cerr << "[" << i << "] " << auxinfo[i].tag[0] << auxinfo[i].tag[1] << std::endl;
					}
				}

				#if 0
				if ( pbam )
				{
					std::cerr << readname << "\t" << libmaus2::bambam::BamAlignmentDecoderBase::getReadName(pbam) << std::endl;
				}
				#endif

				bool const supplementary = (z != il);

				// bool const primary = libmaus2::dazzler::align::OverlapData::getPrimaryFlag(P.first);
				uint64_t const readlen = (package.readsdb)->reads[bread].rlen;

				if ( z==ilow )
				{
					// compute reverse complement
					context.ARC.ensureSize(readlen + 2);

					char * const ra = context.ARC.begin();
					char * rp = ra + readlen + 2;
					char const * src = reinterpret_cast<char const *>((package.readsdb)->bases) + (package.readsdb)->reads[bread].boff - 1;

					*(--rp)	= (*(src++));
					while ( rp != ra+1 )
						*(--rp)	= (*(src++)) ^ 3;
					*(--rp)	= (*(src++));
				}

				uint64_t oauxadd = 0;
				libmaus2::dazzler::align::LASToBamConverterBase::AuxTagIntegerAddRequest req_i("ci",lchainid);
				libmaus2::dazzler::align::LASToBamConverterBase::AuxTagIntegerAddRequest req_n("cn",numchains);
				libmaus2::dazzler::align::LASToBamConverterBase::AuxTagIntegerAddRequest req_j("cj",z-il);
				libmaus2::dazzler::align::LASToBamConverterBase::AuxTagIntegerAddRequest req_l("cl",ih-il);

				Aauxadd.push(oauxadd,&req_i);
				Aauxadd.push(oauxadd,&req_n);
				Aauxadd.push(oauxadd,&req_j);
				Aauxadd.push(oauxadd,&req_l);
				for ( uint64_t i = 0; i < oauxcopy; ++i )
					Aauxadd.push(oauxadd,&Aauxcopy[i]);

				// MD,NM,AS,ci,cn,cj,cl

				libmaus2::dazzler::align::LASToBamConverterBase::AuxTagAddRequest const ** aux_a = &Aauxadd[0];
				libmaus2::dazzler::align::LASToBamConverterBase::AuxTagAddRequest const ** aux_e = &Aauxadd[oauxadd];

				try
				{
					if ( ! inverse )
					{
						context.converter.convert(
							P.first,
							reinterpret_cast<char const *>((package.refdb)->bases) + (package.refdb)->reads[aread].boff,
							(package.refdb)->reads[aread].rlen,
							reinterpret_cast<char const *>((package.readsdb)->bases) + (package.readsdb)->reads[bread].boff,
							(package.readsdb)->reads[bread].rlen,
							readname,
							context.fragment,
							secondary,
							supplementary,
							**(package.Pbamheader),
							aux_a,
							aux_e
						);
					}
					else
					{
						context.converter.convert(
							P.first,
							reinterpret_cast<char const *>((package.refdb)->bases) + (package.refdb)->reads[aread].boff,
							(package.refdb)->reads[aread].rlen,
							context.ARC.begin()+1 /* skip terminator */,
							readlen,
							readname,
							context.fragment,
							secondary,
							supplementary,
							**(package.Pbamheader),
							aux_a,
							aux_e
						);
					}
				}
				catch(std::exception const & ex)
				{
					std::cerr << ex.what() << std::endl;
					std::cerr
						<< aread << " "
						<< bread << " "
						<< flags
						<< std::endl;
				}
			}


			il = ih;
		}

		nextexpt += 1;

		ilow = ihigh;
	}

	package.finsem->post();
}


static void semWait(libmaus2::parallel::PosixSemaphore & finsem, uint64_t const t)
{
	for ( uint64_t i = 0; i < t; ++i )
		finsem.wait();
}

template<typename _package_type>
struct GenericWorkPackage : public libmaus2::parallel::SimpleThreadWorkPackage
{
	typedef _package_type package_type;
	typedef GenericWorkPackage<package_type> this_type;
	typedef typename libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef typename libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	package_type package;
	std::string packagename;

	static std::string getPackageNameStatic()
	{
		return libmaus2::util::Demangle::demangle<this_type>();
	}

	GenericWorkPackage() : libmaus2::parallel::SimpleThreadWorkPackage(), package(), packagename(getPackageNameStatic()) {}
	GenericWorkPackage(
		uint64_t const priority,
		uint64_t const dispatcherid,
		package_type const & rpackage
	) : libmaus2::parallel::SimpleThreadWorkPackage(priority,dispatcherid), package(rpackage), packagename(getPackageNameStatic())
	{}

	char const * getPackageName() const
	{
		return packagename.c_str();
	}
};

struct SearchKmersDispatcher : public libmaus2::parallel::SimpleThreadWorkPackageDispatcher
{
	typedef SearchKmersDispatcher this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	virtual void dispatch(libmaus2::parallel::SimpleThreadWorkPackage * P, libmaus2::parallel::SimpleThreadPoolInterfaceEnqueTermInterface & /* tpi */)
	{
		GenericWorkPackage<PackageSearchKmers> * RP = dynamic_cast< GenericWorkPackage<PackageSearchKmers> * >(P);
		assert ( RP );

		threadSearchKmers(RP->package);
	}
};

struct ConcatSADispatcher : public libmaus2::parallel::SimpleThreadWorkPackageDispatcher
{
	typedef ConcatSADispatcher this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	virtual void dispatch(libmaus2::parallel::SimpleThreadWorkPackage * P, libmaus2::parallel::SimpleThreadPoolInterfaceEnqueTermInterface & /* tpi */)
	{
		GenericWorkPackage<PackageConcatSA> * RP = dynamic_cast< GenericWorkPackage<PackageConcatSA> * >(P);
		assert ( RP );

		threadConcatSA(RP->package);
	}
};

struct FillHistDispatcher : public libmaus2::parallel::SimpleThreadWorkPackageDispatcher
{
	typedef FillHistDispatcher this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	virtual void dispatch(libmaus2::parallel::SimpleThreadWorkPackage * P, libmaus2::parallel::SimpleThreadPoolInterfaceEnqueTermInterface & /* tpi */)
	{
		GenericWorkPackage<PackageFillHist> * RP = dynamic_cast< GenericWorkPackage<PackageFillHist> * >(P);
		assert ( RP );

		threadFillHist(RP->package);
	}
};

struct UnifyPackageSizesDispatcher : public libmaus2::parallel::SimpleThreadWorkPackageDispatcher
{
	typedef UnifyPackageSizesDispatcher this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	virtual void dispatch(libmaus2::parallel::SimpleThreadWorkPackage * P, libmaus2::parallel::SimpleThreadPoolInterfaceEnqueTermInterface & /* tpi */)
	{
		GenericWorkPackage<PackageUnifyPackageSizes> * RP = dynamic_cast< GenericWorkPackage<PackageUnifyPackageSizes> * >(P);
		assert ( RP );

		threadUnifyPackageSizes(RP->package);
	}
};

struct MergeHistogramsDispatcher : public libmaus2::parallel::SimpleThreadWorkPackageDispatcher
{
	typedef MergeHistogramsDispatcher this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	virtual void dispatch(libmaus2::parallel::SimpleThreadWorkPackage * P, libmaus2::parallel::SimpleThreadPoolInterfaceEnqueTermInterface & /* tpi */)
	{
		GenericWorkPackage<PackageMergeHistograms> * RP = dynamic_cast< GenericWorkPackage<PackageMergeHistograms> * >(P);
		assert ( RP );

		threadMergeHistograms(RP->package);
	}
};

struct ConcatenateUniqueDispatcher : public libmaus2::parallel::SimpleThreadWorkPackageDispatcher
{
	typedef ConcatenateUniqueDispatcher this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	virtual void dispatch(libmaus2::parallel::SimpleThreadWorkPackage * P, libmaus2::parallel::SimpleThreadPoolInterfaceEnqueTermInterface & /* tpi */)
	{
		GenericWorkPackage<PackageConcatenateUnique> * RP = dynamic_cast< GenericWorkPackage<PackageConcatenateUnique> * >(P);
		assert ( RP );

		threadConcatenateUnique(RP->package);
	}
};

struct SumIntervalsDispatcher : public libmaus2::parallel::SimpleThreadWorkPackageDispatcher
{
	typedef SumIntervalsDispatcher this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	virtual void dispatch(libmaus2::parallel::SimpleThreadWorkPackage * P, libmaus2::parallel::SimpleThreadPoolInterfaceEnqueTermInterface & /* tpi */)
	{
		GenericWorkPackage<PackageSumIntervals> * RP = dynamic_cast< GenericWorkPackage<PackageSumIntervals> * >(P);
		assert ( RP );
		threadSumIntervals(RP->package);
	}
};

struct PrefixSumsDispatcher : public libmaus2::parallel::SimpleThreadWorkPackageDispatcher
{
	typedef PrefixSumsDispatcher this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	virtual void dispatch(libmaus2::parallel::SimpleThreadWorkPackage * P, libmaus2::parallel::SimpleThreadPoolInterfaceEnqueTermInterface & /* tpi */)
	{
		GenericWorkPackage<PackagePrefixSums> * RP = dynamic_cast< GenericWorkPackage<PackagePrefixSums> * >(P);
		assert ( RP );
		threadPrefixSums(RP->package);
	}
};

struct LookupSADispatcher : public libmaus2::parallel::SimpleThreadWorkPackageDispatcher
{
	typedef LookupSADispatcher this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	virtual void dispatch(libmaus2::parallel::SimpleThreadWorkPackage * P, libmaus2::parallel::SimpleThreadPoolInterfaceEnqueTermInterface & /* tpi */)
	{
		GenericWorkPackage<PackageLookupSA> * RP = dynamic_cast< GenericWorkPackage<PackageLookupSA> * >(P);
		assert ( RP );
		threadLookupSA(RP->package);
	}
};

struct CompactPosDispatcher : public libmaus2::parallel::SimpleThreadWorkPackageDispatcher
{
	typedef CompactPosDispatcher this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	virtual void dispatch(libmaus2::parallel::SimpleThreadWorkPackage * P, libmaus2::parallel::SimpleThreadPoolInterfaceEnqueTermInterface & /* tpi */)
	{
		GenericWorkPackage<PackageCompactPos> * RP = dynamic_cast< GenericWorkPackage<PackageCompactPos> * >(P);
		assert ( RP );
		threadCompactPos(RP->package);
	}
};

struct ConcatQueryKmersDispatcher : public libmaus2::parallel::SimpleThreadWorkPackageDispatcher
{
	typedef ConcatQueryKmersDispatcher this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	virtual void dispatch(libmaus2::parallel::SimpleThreadWorkPackage * P, libmaus2::parallel::SimpleThreadPoolInterfaceEnqueTermInterface & /* tpi */)
	{
		GenericWorkPackage<PackageConcatQueryKmers> * RP = dynamic_cast< GenericWorkPackage<PackageConcatQueryKmers> * >(P);
		assert ( RP );
		threadConcatQueryKmers(RP->package);
	}
};

struct EncodeBamDispatcher : public libmaus2::parallel::SimpleThreadWorkPackageDispatcher
{
	typedef EncodeBamDispatcher this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	virtual void dispatch(libmaus2::parallel::SimpleThreadWorkPackage * P, libmaus2::parallel::SimpleThreadPoolInterfaceEnqueTermInterface & /* tpi */)
	{
		GenericWorkPackage<PackageEncodeBam> * RP = dynamic_cast< GenericWorkPackage<PackageEncodeBam> * >(P);
		assert ( RP );
		threadEncodeBam(RP->package);
	}
};

struct CompressBgzfDispatcher : public libmaus2::parallel::SimpleThreadWorkPackageDispatcher
{
	typedef EncodeBamDispatcher this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	virtual void dispatch(libmaus2::parallel::SimpleThreadWorkPackage * P, libmaus2::parallel::SimpleThreadPoolInterfaceEnqueTermInterface & /* tpi */)
	{
		GenericWorkPackage<PackageCompressBgzf> * RP = dynamic_cast< GenericWorkPackage<PackageCompressBgzf> * >(P);
		assert ( RP );
		RP->package.dispatch();
	}
};

#if 0
enum dispatcher_ids
{
	SearchKmersDispatcher_id = 0,
	ConcatSADispatcher_id = 1,
	FillHistDispatcher_id = 2,
	UnifyPackageSizesDispatcher_id = 3,
	MergeHistogramsDispatcher_id = 4,
	ConcatenateUniqueDispatcher_id = 5,
	SumIntervalsDispatcher_id = 6,
	PrefixSumsDispatcher_id = 7,
	LookupSADispatcher_id = 8,
	CompactPosDispatcher_id = 9,
	ConcatQueryKmersDispatcher_id = 10,
	EncodeBamDispatcher_id = 11
};
#endif

int damapper_bwt(libmaus2::util::ArgParser const & arg)
{
	// default k
	unsigned int const defk = 20;
	unsigned int const k = arg.argPresent("k") ? arg.getUnsignedNumericArg<uint64_t>("k") : defk;

	// verbosity
	unsigned int const defv = 1;
	unsigned int const verbose = arg.argPresent("v") ? arg.getUnsignedNumericArg<uint64_t>("v") : defv;

	// maximum number of input bases per block
	uint64_t const defo = 256*1024*1024;
	uint64_t const maxo = arg.argPresent("i") ? std::max(static_cast<uint64_t>(1),arg.getUnsignedNumericArg<uint64_t>("i")) : defo;

	// length for K-mer cache
	unsigned int const defcachek = 12;
	unsigned int const cache_k = arg.uniqueArgPresent("K") ? arg.getUnsignedNumericArg<uint64_t>("K") : defcachek;

	MEM_PHYSICAL = libmaus2::util::MemoryStatistics::getPhysicalMemory();
	MEM_LIMIT = MEM_PHYSICAL;

	uint64_t const defconstrsasamplingrate = 4;
	uint64_t const constrsasamplingrate = arg.argPresent("sasamplingrate") ? std::max(static_cast<uint64_t>(1),arg.getUnsignedNumericArg<uint64_t>("sasamplingrate")) : defconstrsasamplingrate;
	uint64_t const deficonstrsasamplingrate = 32;
	uint64_t const iconstrsasamplingrate = arg.argPresent("isasamplingrate") ? std::max(static_cast<uint64_t>(1),arg.getUnsignedNumericArg<uint64_t>("isasamplingrate")) : deficonstrsasamplingrate;
	uint64_t const defaultbwtconstrmem = (MEM_PHYSICAL*3)/4;
	uint64_t const bwtconstrmem = arg.argPresent("bwtconstrmem") ? std::max(static_cast<uint64_t>(1),arg.getUnsignedNumericArg<uint64_t>("bwtconstrmem")) : defaultbwtconstrmem;

	int const level = arg.uniqueArgPresent("z") ? arg.getParsedArg<int64_t>("z") : getDefaultLevel();
	BgzfDeflateOutputBufferBaseAllocator bgzfalloc(level);
	libmaus2::parallel::LockedGrowingFreeList<libmaus2::lz::BgzfDeflateOutputBufferBase,BgzfDeflateOutputBufferBaseAllocator,BgzfDeflateOutputBufferBaseTypeInfo> bgzffreelist(bgzfalloc);
	BgzfDeflateZStreamBaseAllocator bgzfzalloc(level);
	libmaus2::parallel::LockedGrowingFreeList<libmaus2::lz::BgzfDeflateZStreamBase,BgzfDeflateZStreamBaseAllocator,BgzfDeflateZStreamBaseTypeInfo> zstreambasefreelist(bgzfzalloc);

	#if defined(_OPENMP)
	uint64_t const defnumproc = libmaus2::parallel::NumCpus::getNumLogicalProcessors();
	uint64_t const numthreads = arg.argPresent("p") ? std::max(static_cast<uint64_t>(1),arg.getUnsignedNumericArg<uint64_t>("p")) : defnumproc;
	#else
	uint64_t const numthreads = 1;
	#endif

	libmaus2::parallel::SimpleThreadPool STP(numthreads);

	try
	{
		SearchKmersDispatcher SKD;
		ConcatSADispatcher CSAD;
		FillHistDispatcher FHD;
		UnifyPackageSizesDispatcher UPSD;
		MergeHistogramsDispatcher MHD;
		ConcatenateUniqueDispatcher CUD;
		SumIntervalsDispatcher SID;
		PrefixSumsDispatcher PSD;
		LookupSADispatcher LSAD;
		CompactPosDispatcher CPD;
		ConcatQueryKmersDispatcher CQKD;
		EncodeBamDispatcher EBD;
		CompressBgzfDispatcher CBD;
		uint64_t const SearchKmersDispatcher_id = STP.getNextDispatcherId();
		uint64_t const ConcatSADispatcher_id = STP.getNextDispatcherId();
		uint64_t const FillHistDispatcher_id = STP.getNextDispatcherId();
		uint64_t const UnifyPackageSizesDispatcher_id = STP.getNextDispatcherId();
		uint64_t const MergeHistogramsDispatcher_id = STP.getNextDispatcherId();
		uint64_t const ConcatenateUniqueDispatcher_id = STP.getNextDispatcherId();
		uint64_t const SumIntervalsDispatcher_id = STP.getNextDispatcherId();
		uint64_t const PrefixSumsDispatcher_id = STP.getNextDispatcherId();
		uint64_t const LookupSADispatcher_id = STP.getNextDispatcherId();
		uint64_t const CompactPosDispatcher_id = STP.getNextDispatcherId();
		uint64_t const ConcatQueryKmersDispatcher_id = STP.getNextDispatcherId();
		uint64_t const EncodeBamDispatcher_id = STP.getNextDispatcherId();
		uint64_t const CompressBgzfDispatcher_id = STP.getNextDispatcherId();
		STP.registerDispatcher(SearchKmersDispatcher_id,&SKD);
		STP.registerDispatcher(ConcatSADispatcher_id,&CSAD);
		STP.registerDispatcher(FillHistDispatcher_id,&FHD);
		STP.registerDispatcher(UnifyPackageSizesDispatcher_id,&UPSD);
		STP.registerDispatcher(MergeHistogramsDispatcher_id,&MHD);
		STP.registerDispatcher(ConcatenateUniqueDispatcher_id,&CUD);
		STP.registerDispatcher(SumIntervalsDispatcher_id,&SID);
		STP.registerDispatcher(PrefixSumsDispatcher_id,&PSD);
		STP.registerDispatcher(LookupSADispatcher_id,&LSAD);
		STP.registerDispatcher(CompactPosDispatcher_id,&CPD);
		STP.registerDispatcher(ConcatQueryKmersDispatcher_id,&CQKD);
		STP.registerDispatcher(EncodeBamDispatcher_id,&EBD);
		STP.registerDispatcher(CompressBgzfDispatcher_id,&CBD);

		std::string const s_supstrat = arg.uniqueArgPresent("S") ? arg["S"] : "none";

		libmaus2::dazzler::align::LASToBamConverterBase::supplementary_seq_strategy_t e_supstrat;

		if ( s_supstrat == "none" )
		{
			e_supstrat = libmaus2::dazzler::align::LASToBamConverterBase::supplementary_seq_strategy_none;
		}
		else if ( s_supstrat == "soft" )
		{
			e_supstrat = libmaus2::dazzler::align::LASToBamConverterBase::supplementary_seq_strategy_soft;
		}
		else if ( s_supstrat == "hard" )
		{
			e_supstrat = libmaus2::dazzler::align::LASToBamConverterBase::supplementary_seq_strategy_hard;
		}
		else
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "[E] unknown storage strategy " << s_supstrat << std::endl;
			lme.finish();
			throw lme;
		}

		// std::cerr << "e_supstrat=" << static_cast<int>(e_supstrat) << " " << s_supstrat << std::endl;

		uint64_t const defmemlimit = 0;
		uint64_t const memlimit = arg.uniqueArgPresent("M") ? arg.getUnsignedNumericArg<uint64_t>("M") : defmemlimit;
		if ( memlimit )
			MEM_LIMIT = memlimit;

		if ( verbose > 1 )
		{
			std::cerr << "[V] physical memory " << MEM_PHYSICAL << std::endl;
			if ( memlimit )
				std::cerr << "[V] mem limit " << memlimit << std::endl;
		}

		std::string const fastaname = arg[0];
		std::string const combfn = arg.uniqueArgPresent("Q") ? arg["Q"] : (fastaname + ".damapper_bwt");
		std::string const readsformat = arg.uniqueArgPresent("I") ? arg["I"] : "fasta";

		if (
			(! libmaus2::util::GetFileSize::fileExists(combfn))
			||
			libmaus2::util::GetFileSize::isOlder(combfn,fastaname)
		)
		{
			DNAIndexBuild forwindex;
			DNAIndexBuild rcindex;

			std::string const tmpprefix = arg.uniqueArgPresent("T") ? arg["T"] : libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname);
			std::string const indexcreatefn = tmpprefix + "_damapper_bwt_index_create";
			libmaus2::util::TempFileRemovalContainer::addTempFile(indexcreatefn);

			libmaus2::aio::OutputStreamInstance::unique_ptr_type combOSI(new libmaus2::aio::OutputStreamInstance(indexcreatefn));
			forwindex.setupFromFasta(fastaname,tmpprefix,verbose,constrsasamplingrate,iconstrsasamplingrate,bwtconstrmem,numthreads,*combOSI);
			rcindex.setupFromForwardCompact(fastaname,tmpprefix,forwindex.fainame,forwindex.compactname,forwindex.metaname,forwindex.replname,verbose,constrsasamplingrate,iconstrsasamplingrate,bwtconstrmem,numthreads,*combOSI);

			std::string const m5info = tmpprefix+".m5infos";
			libmaus2::aio::InputStreamInstance ISI(fastaname);
			libmaus2::fastx::FastAStreamSet FASS(ISI);
			std::map<std::string,std::string> m5map = FASS.computeMD5(false,false);

			libmaus2::util::NumberSerialisation::serialiseNumber(*combOSI,m5map.size());
			for ( std::map<std::string,std::string>::const_iterator ita = m5map.begin(); ita != m5map.end(); ++ita )
			{
				libmaus2::util::StringSerialisation::serialiseString(*combOSI,ita->first);
				libmaus2::util::StringSerialisation::serialiseString(*combOSI,ita->second);
			}

			combOSI->flush();
			combOSI.reset();

			forwindex.cleanup();
			rcindex.cleanup();

			// rename
			libmaus2::aio::OutputStreamFactoryContainer::rename(indexcreatefn,combfn);
		}

		DNAIndex forwindex;
		DNAIndex rcindex;

		libmaus2::aio::InputStreamInstance::unique_ptr_type combISI(new libmaus2::aio::InputStreamInstance(combfn));
		forwindex.loadFromSingleForward(*combISI,numthreads,cache_k,verbose);
		rcindex.loadFromSingleReverseComplement(*combISI,numthreads,cache_k,verbose);

		std::map<std::string,std::string> m5map;
		uint64_t const sm5 = libmaus2::util::NumberSerialisation::deserialiseNumber(*combISI);
		for ( uint64_t i = 0; i < sm5; ++i )
		{
			std::string const key = libmaus2::util::StringSerialisation::deserialiseString(*combISI);
			std::string const value = libmaus2::util::StringSerialisation::deserialiseString(*combISI);
			m5map[key] = value;
		}

		combISI.reset();


		std::vector < DNAIndex * > indexes;
		indexes.push_back(&forwindex);
		indexes.push_back(&rcindex);

		float freq[4] = { 0.25, 0.25, 0.25, 0.25 };
		int64_t const tspace = 100;
		Align_Spec * aspec = New_Align_Spec(0.85, tspace, &freq[0]);
		Set_Filter_Params(k, 0, numthreads);

		libmaus2::fastx::LineBufferFastAReader::unique_ptr_type LBFA;
		libmaus2::bambam::BamAlignmentDecoderWrapper::unique_ptr_type BAD;

		if ( readsformat == "fasta" )
		{
			libmaus2::fastx::LineBufferFastAReader::unique_ptr_type TLBFA(new libmaus2::fastx::LineBufferFastAReader(std::cin));
			LBFA = UNIQUE_PTR_MOVE(TLBFA);
		}
		else if ( readsformat == "bam" )
		{
			libmaus2::bambam::BamAlignmentDecoderWrapper::unique_ptr_type TBAD(libmaus2::bambam::BamAlignmentDecoderFactory::construct());
			BAD = UNIQUE_PTR_MOVE(TBAD);
		}
		else
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "[E] unsupported input format " << readsformat << std::endl;
			lme.finish();
			throw lme;
		}

		libmaus2::bambam::BamHeader const * inheader = BAD ? (&(BAD->getDecoder().getHeader())) : 0;

		bool running = true;

		libmaus2::timing::RealTimeClock roundclock;
		libmaus2::timing::RealTimeClock lclock;

		// read names
		libmaus2::autoarray::AutoArray<char> Areadnames;
		// read base data
		libmaus2::autoarray::AutoArray<char> Areaddata;
		// BAM data (if any)
		libmaus2::autoarray::AutoArray<uint8_t> Abam;
		// read meta data
		libmaus2::autoarray::AutoArray< libmaus2::fastx::LineBufferFastAReader::ReadMeta > O;
		// bam meta data
		libmaus2::autoarray::AutoArray< uint64_t > Obam;

		uint64_t const refnumposbytes = forwindex.getNumPosBytes();
		uint64_t const refnumseqbytes = forwindex.getNumSeqBytes();
		uint64_t const numkeybits = 2*k;
		uint64_t const numkeybytes = (numkeybits + 7)/8;

		std::vector<unsigned int> qfullkeybytes;
		#if defined(LIBMAUS2_BYTE_ORDER_LITTLE_ENDIAN)
		// position in sequence
		for ( unsigned int i = 0; i < refnumposbytes; ++i )
			qfullkeybytes.push_back(sizeof(uint64_t)+i);
		// sequence id
		for ( unsigned int i = 0; i < refnumseqbytes; ++i )
			qfullkeybytes.push_back(sizeof(uint64_t)+sizeof(uint32_t)+i);
		// code
		for ( unsigned int i = 0; i < numkeybytes; ++i )
			qfullkeybytes.push_back(i);
		#else
		#error "Unsupported byte order"
		#endif

		uint64_t const packsperthread = 32;
		uint64_t const tnumpacks = packsperthread*numthreads;

		// query intervals on suffix array (divided into packages)
		libmaus2::autoarray::AutoArray< libmaus2::autoarray::AutoArray < RefKMer >::unique_ptr_type > RKM(tnumpacks);
		libmaus2::autoarray::AutoArray< libmaus2::autoarray::AutoArray < QueryKMer >::unique_ptr_type > QKM(tnumpacks);

		libmaus2::autoarray::AutoArray< uint64_t > NUMRK(tnumpacks+1);
		libmaus2::autoarray::AutoArray< uint64_t > NUMQK(tnumpacks+1);
		libmaus2::autoarray::AutoArray< uint64_t > NUMRV(tnumpacks+1);

		// ref db read array
		libmaus2::autoarray::AutoArray< HITS_READ > AREFHITREADS;
		libmaus2::autoarray::AutoArray< HITS_READ > AREADREADS;
		HITS_DB refdb;

		for ( uint64_t i = 0; i < RKM.size(); ++i )
		{
			libmaus2::autoarray::AutoArray < RefKMer >::unique_ptr_type tptr(new libmaus2::autoarray::AutoArray < RefKMer >);
			RKM[i] = UNIQUE_PTR_MOVE(tptr);
		}
		for ( uint64_t i = 0; i < QKM.size(); ++i )
		{
			libmaus2::autoarray::AutoArray < QueryKMer >::unique_ptr_type tptr(new libmaus2::autoarray::AutoArray < QueryKMer >);
			QKM[i] = UNIQUE_PTR_MOVE(tptr);
		}

		libmaus2::autoarray::AutoArray < RefKMer > GRKM;
		libmaus2::autoarray::AutoArray < RefKMer > TGRKM;
		libmaus2::autoarray::AutoArray < QueryKMer, libmaus2::autoarray::alloc_type_c > GQKM;
		libmaus2::autoarray::AutoArray < QueryKMer, libmaus2::autoarray::alloc_type_c > TGQKM;
		libmaus2::autoarray::AutoArray < QueryKMer, libmaus2::autoarray::alloc_type_c > GRQKM;
		libmaus2::autoarray::AutoArray < QueryKMer, libmaus2::autoarray::alloc_type_c > TGRQKM;

		// key bytes for code sorting
		std::vector<unsigned int> keybytes;
		#if defined(LIBMAUS2_BYTE_ORDER_LITTLE_ENDIAN)
		for ( uint64_t i = 0; i < numkeybytes; ++i )
			keybytes.push_back(i);
		#elif defined(LIBMAUS2_BYTE_ORDER_BIG_ENDIAN)
		for ( uint64_t i = 0; i < numkeybytes; ++i )
			keybytes.push_back(sizeof(uint64_t) - i - 1);
		#else
		#error "Unknown byte order"
		#endif

		// invalid query kmer
		QueryKMer QIff;
		QIff.code = std::numeric_limits<uint64_t>::max();
		QIff.rpos = std::numeric_limits<uint32_t>::max();
		QIff.readid = std::numeric_limits<uint32_t>::max();
		QueryKMer QI0;
		QI0.code = 0;
		QI0.rpos = std::numeric_limits<uint32_t>::max();
		QI0.readid = std::numeric_limits<uint32_t>::max();

		libmaus2::bambam::BamHeader::unique_ptr_type Pbamheader(forwindex.getBamHeader(m5map,fastaname,inheader,arg));
		libmaus2::lz::BgzfDeflate<std::ostream>::unique_ptr_type bgzfout(new libmaus2::lz::BgzfDeflate<std::ostream>(std::cout,level));
		Pbamheader->serialise(*bgzfout);
		bgzfout->flush();

		AlignerAllocator alignalloc;
		libmaus2::autoarray::AutoArray< libmaus2::lcs::Aligner::shared_ptr_type > Aaligners(numthreads);
		for ( uint64_t i = 0; i < Aaligners.size(); ++i )
		{
			libmaus2::lcs::Aligner::shared_ptr_type Paligner = alignalloc();
			Aaligners[i] = Paligner;
		}

		libmaus2::autoarray::AutoArray< libmaus2::lcs::AlignmentTraceContainer::unique_ptr_type > AATC(numthreads);
		for ( uint64_t i = 0; i < AATC.size(); ++i )
		{
			libmaus2::lcs::AlignmentTraceContainer::unique_ptr_type Tptr(new libmaus2::lcs::AlignmentTraceContainer);
			AATC[i] = UNIQUE_PTR_MOVE(Tptr);
		}


		libmaus2::autoarray::AutoArray<LASToBamContext::unique_ptr_type> Acontexts(numthreads);
		libmaus2::dazzler::align::RefMapEntryVector refmap;
		for ( uint64_t i = 0; i < Pbamheader->getNumRef(); ++i )
			refmap.push_back(libmaus2::dazzler::align::RefMapEntry(i,0));
		for ( uint64_t i = 0; i < numthreads; ++i )
		{
			LASToBamContext::unique_ptr_type Tptr(new LASToBamContext(
				tspace,true /* mdnm */,
				e_supstrat,
				std::string(),refmap));
			Acontexts[i] = UNIQUE_PTR_MOVE(Tptr);
		}

		uint64_t const histthres = 64*1024;
		uint64_t const histmult = histthres + 1;
		libmaus2::autoarray::AutoArray< libmaus2::autoarray::AutoArray<uint64_t>::unique_ptr_type> Ahist(numthreads);
		for ( uint64_t i = 0; i < numthreads; ++i )
		{
			libmaus2::autoarray::AutoArray<uint64_t>::unique_ptr_type thist(new libmaus2::autoarray::AutoArray<uint64_t>(histmult));
			Ahist[i] = UNIQUE_PTR_MOVE(thist);
		}

		libmaus2::timing::RealTimeClock accclock;
		accclock.start();

		uint64_t totalreads = 0;

		libmaus2::parallel::PosixSemaphore finsem;

		for ( uint64_t readpart = 1 ; running; ++readpart )
		{
			roundclock.start();

			uint64_t o_name = 0;
			uint64_t o_data = 0;
			uint64_t prev_o_name = o_name;
			uint64_t rlen = 0;
			bool const termdata = false;
			bool const mapdata = true;
			bool const paddata = true;
			char const padsym = 4;
			bool havedata = true;
			uint64_t oo = 0;
			uint64_t oobam = 0;

			lclock.start();
			if ( LBFA )
			{
				// load read block
				while ( o_data < maxo && (havedata = LBFA->getNext(Areadnames,o_name,Areaddata,o_data,rlen,termdata,mapdata,paddata,padsym)) )
				{
					#if 0
					std::cerr << data.begin()+o_name << std::endl;
					std::cerr.write(
						data.begin() + o_data,
						rlen
					);
					std::cerr.put('\n');
					#endif

					uint64_t readoffset = o_data - rlen - 2*(paddata ? 1 : 0) + (paddata ? 1 : 0);

					assert ( Areaddata [ readoffset - 1 ] == padsym );
					assert ( Areaddata [ readoffset + rlen ] == padsym );

					O.push(oo,libmaus2::fastx::LineBufferFastAReader::ReadMeta(prev_o_name,readoffset,rlen));
					prev_o_name = o_name;
				}
			}
			else if ( BAD )
			{
				libmaus2::bambam::BamAlignmentDecoder & decoder = BAD->getDecoder();
				libmaus2::bambam::BamAlignment const & algn = decoder.getAlignment();
				uint64_t o_bam = 0;

				while ( o_data < maxo && (havedata=decoder.readAlignment()) )
				{
					char const * cname = algn.getName();
					uint64_t const lname = algn.getLReadName();
					uint64_t const rlen = algn.getLseq();
					assert ( strlen(cname)+1 == lname );

					Areadnames.ensureSize(o_name + lname);
					std::copy(cname,cname+lname,Areadnames.begin() + o_name);

					Areaddata.ensureSize(o_data + 2 + rlen);
					Areaddata [ o_data ] = padsym;
					Areaddata [ o_data + rlen + 1 ] = padsym;
					algn.decodeRead(Areaddata.begin() + o_data + 1, rlen);

					for ( uint64_t i = 0; i < rlen; ++i )
						Areaddata [ o_data + 1 + i ] = libmaus2::fastx::mapChar(Areaddata [ o_data + 1 + i ]);

					uint64_t const readoffset = o_data + 1;

					assert ( Areaddata [ readoffset - 1 ] == padsym );
					assert ( Areaddata [ readoffset + rlen ] == padsym );

					Abam.ensureSize(o_bam + algn.blocksize);
					std::copy(
						algn.D.begin(),
						algn.D.begin()+algn.blocksize,
						Abam.begin() + o_bam
					);

					O.push(oo,libmaus2::fastx::LineBufferFastAReader::ReadMeta(o_name,readoffset,rlen));
					Obam.push(oobam,o_bam);

					o_name += lname;
					o_data += rlen + 2;
					o_bam += algn.blocksize;
				}
				Obam.push(oobam,o_bam);
			}

			uint64_t const numreads = oo;
			if ( verbose > 1 )
				std::cerr << "[V] got " << numreads << " reads in time " << lclock.formatTime(lclock.getElapsedSeconds()) << std::endl;

			std::vector<std::string> readnames;
			for ( uint64_t i = 0; i < numreads; ++i )
				readnames.push_back(Areadnames.begin() + O[i].nameoff);

			// set up HITS_DB datastructure for reads
			uint64_t readmaxlen = 0;
			uint64_t readtotlen = 0;
			for ( uint64_t i = 0; i < numreads; ++i )
			{
				readmaxlen = std::max(readmaxlen,O[i].len);
				readtotlen += O[i].len;
			}

			#if 0
			uint64_t const numreadlenbits = readmaxlen ? libmaus2::math::numbits(readmaxlen-1) : 0;
			uint64_t const numreadlenbytes = (numreadlenbits + 7)/8;
			#endif

			AREADREADS.ensureSize(numreads + 2);
			static char readpath[] = "reads.db";

			HITS_DB readsdb;
			readsdb.ureads = numreads;
			readsdb.treads = numreads;
			readsdb.cutoff = -1;
			// readsdb.all = 0;
			readsdb.freq[0] = readsdb.freq[1] = readsdb.freq[2] = readsdb.freq[3] = 0.25;
			readsdb.maxlen = readmaxlen;
			readsdb.totlen = readtotlen;
			readsdb.nreads = numreads;
			readsdb.trimmed = 1;
			readsdb.part = readpart;
			readsdb.ufirst = 0;
			readsdb.tfirst = 0;
			readsdb.path = &readpath[0];
			readsdb.loaded = 1;
			readsdb.bases = Areaddata.begin();
			readsdb.reads = AREADREADS.begin() + 1;
			readsdb.tracks = NULL;

			for ( uint64_t i = 0; i < numreads; ++i )
			{
				readsdb.reads[i].origin = i;
				readsdb.reads[i].rlen = O[i].len;
				readsdb.reads[i].fpulse = 0;
				readsdb.reads[i].boff = O[i].dataoff;
				readsdb.reads[i].coff = -1;
				readsdb.reads[i].flags = DB_BEST;
			}

			(reinterpret_cast<int*>(readsdb.reads))[-1] = numreads;
			(reinterpret_cast<int*>(readsdb.reads))[-2] = numreads;

			// number of non empty intervals per read
			#if 0
			RI.ensureSize(numreads);
			RS.ensureSize(numreads+1);
			RSN.ensureSize(numreads+1);
			RSE.ensureSize(numreads+1);
			#endif

			uint64_t const readsperpack = (numreads + tnumpacks - 1)/tnumpacks;
			uint64_t const numpacks = readsperpack ? (numreads + readsperpack - 1)/readsperpack : 0;
			assert ( numpacks <= tnumpacks );

			struct KmerPos
			{
				// kmer code
				uint64_t code;
				// position in read
				int32_t rpos;
				// read id
				int32_t read;
			};

			uint64_t indexyield = 0;

			for ( uint64_t index_i = 0; index_i < indexes.size(); ++index_i )
			{
				DNAIndex & index = *indexes[index_i];

				// #define IDEBUG
				#if defined(IDEBUG)
				libmaus2::bitio::CompactArray::unique_ptr_type CTEXT(index.loadText());
				#endif

				// search kmers
				lclock.start();
				std::vector < GenericWorkPackage<PackageSearchKmers> > VSAR(numpacks);
				for ( uint64_t t = 0; t < numpacks; ++t )
				{
					VSAR[t] = GenericWorkPackage<PackageSearchKmers>(0,SearchKmersDispatcher_id,PackageSearchKmers(t,readsperpack,numreads,k,&RKM,&QKM,&Areaddata,&O,&NUMRK,&NUMQK,&index,&finsem));
					STP.enque(&VSAR[t]);
				}
				semWait(finsem,numpacks);
				if ( verbose > 1 )
					std::cerr << "[V] computed SA ranges in time " << lclock.formatTime(lclock.getElapsedSeconds()) << std::endl;

				// concatenate reference kmer sequences
				lclock.start();
				libmaus2::util::PrefixSums::prefixSums(NUMRK.begin(),NUMRK.begin()+numpacks+1);
				GRKM.ensureSize(NUMRK[numpacks]);

				std::vector < GenericWorkPackage<PackageConcatSA> > VCSAD(numpacks);
				for ( uint64_t t = 0; t < numpacks; ++t )
				{
					VCSAD[t] = GenericWorkPackage<PackageConcatSA>(0,ConcatSADispatcher_id,PackageConcatSA(t,&GRKM,&NUMRK,&RKM,&finsem));
					STP.enque(&VCSAD[t]);
				}
				semWait(finsem,numpacks);
				if ( verbose > 1 )
					std::cerr << "[V] concatenated SA ranges in time " << lclock.formatTime(lclock.getElapsedSeconds()) << std::endl;

				// sort suffix ranges by code field (so we can unify the ranges)
				lclock.start();
				TGRKM.ensureSize(NUMRK[numpacks]);
				libmaus2::sorting::InterleavedRadixSort::byteradixsortKeyBytes(
					GRKM.begin(),GRKM.begin()+NUMRK[numpacks],
					TGRKM.begin(),TGRKM.begin()+NUMRK[numpacks],
					numthreads,
					keybytes.begin(),numkeybytes,STP);
				if ( verbose > 1 )
					std::cerr << "[V] sorted concatenated SA ranges in time " << lclock.formatTime(lclock.getElapsedSeconds()) << std::endl;
				if ( numkeybytes % 2 == 1 )
					std::swap(GRKM,TGRKM);

				#if 0
				for ( uint64_t i = 1; i < NUMRK[numpacks]; ++i )
					assert ( GRKM[i-1].code <= GRKM[i].code );
				#endif

				lclock.start();
				// look for code bounds
				std::vector<uint64_t> ubounds;
				uint64_t const upacksize = (NUMRK[numpacks] + numthreads - 1)/numthreads;
				uint64_t ulow = 0;
				while ( ulow < NUMRK[numpacks] )
				{
					ubounds.push_back(ulow);
					uint64_t uhigh = std::min(NUMRK[numpacks],ulow+upacksize);
					while ( uhigh+1 < NUMRK[numpacks] && GRKM[uhigh].code == GRKM[uhigh-1].code )
						++uhigh;

					ulow = uhigh;
				}
				ubounds.push_back(ulow);
				assert ( ulow == NUMRK[numpacks] );

				// check bounds
				for ( uint64_t i = 0; i < ubounds.size(); ++i )
					assert ( i == 0 || (i+1)==ubounds.size() || GRKM[ubounds[i]-1].code != GRKM[ubounds[i]].code );

				std::vector < GenericWorkPackage<PackageFillHist> > VPFH(numthreads);
				for ( uint64_t t = 0; t < numthreads; ++t )
				{
					VPFH[t] = GenericWorkPackage<PackageFillHist>(0,FillHistDispatcher_id,PackageFillHist(t,&Ahist,&finsem));
					STP.enque(&VPFH[t]);
				}
				semWait(finsem,numthreads);

				// unified package sizes
				std::vector<uint64_t> usize(ubounds.size());
				uint64_t const numuni = ubounds.size() ? (ubounds.size()-1) : 0;
				std::vector < GenericWorkPackage<PackageUnifyPackageSizes> > VUPS(numuni);
				for ( uint64_t t = 1; t < ubounds.size(); ++t )
				{
					VUPS[t-1] = GenericWorkPackage<PackageUnifyPackageSizes>(0,UnifyPackageSizesDispatcher_id,PackageUnifyPackageSizes(t,&Ahist,&GRKM,&ubounds,&usize,histthres,&finsem));
					STP.enque(&VUPS[t-1]);
				}
				semWait(finsem,numuni);
				// compute prefix sums
				libmaus2::util::PrefixSums::prefixSums(usize.begin(),usize.end());

				// merge histograms
				uint64_t const tperthread = (histmult + numthreads - 1)/numthreads;
				std::vector < GenericWorkPackage < PackageMergeHistograms > > VPMH(numthreads);
				for ( uint64_t t = 0; t < numthreads; ++t )
				{
					VPMH[t] = GenericWorkPackage < PackageMergeHistograms >(0,MergeHistogramsDispatcher_id,PackageMergeHistograms(t,tperthread,histmult,numthreads,&Ahist,&finsem));
					STP.enque(&VPMH[t]);
				}
				semWait(finsem,numthreads);

				#if 0
				for ( uint64_t i = 0; i < histmult; ++i )
					if ( (*(Ahist[0]))[i] )
					{
						std::cerr << "[H] " << i << " " << (*(Ahist[0]))[i] << std::endl;
					}
				#endif

				// concatenate unique
				uint64_t volatile gc = 0;
				libmaus2::parallel::PosixSpinLock gclock;
				uint64_t const numus = usize.size() ? (usize.size()-1) : 0;
				std::vector < GenericWorkPackage < PackageConcatenateUnique > > VPCU(numus);
				for ( uint64_t t = 1; t < usize.size(); ++t )
				{
					VPCU[t-1] = GenericWorkPackage < PackageConcatenateUnique >(0,ConcatenateUniqueDispatcher_id,PackageConcatenateUnique(t,&GRKM,&TGRKM,&ubounds,&usize,&gclock,&gc,&finsem));
					STP.enque(&VPCU[t-1]);
				}
				semWait(finsem,numus);
				GRKM.swap(TGRKM);
				if ( verbose > 1 )
					std::cerr << "[V] uniquified in time " << lclock.formatTime(lclock.getElapsedSeconds()) << " num out " << usize.back() << " sum " << gc << std::endl;

				lclock.start();
				libmaus2::autoarray::AutoArray<uint64_t,libmaus2::autoarray::alloc_type_c> PP(usize.back()+1,false);
				uint64_t const supacksize = (usize.back() + numthreads - 1)/numthreads;
				uint64_t const suthreads = supacksize ? ((usize.back() + supacksize - 1) /supacksize) : 0;
				PP[0] = 0;

				std::vector < GenericWorkPackage < PackageSumIntervals > > VSI(suthreads);
				for ( uint64_t t = 0; t < suthreads; ++t )
				{
					VSI[t] = GenericWorkPackage < PackageSumIntervals >(0,SumIntervalsDispatcher_id,PackageSumIntervals(t,supacksize,&usize,&GRKM,&PP,&finsem));
					STP.enque(&VSI[t]);
				}
				semWait(finsem,suthreads);

				for ( uint64_t t = 0; t < suthreads; ++t )
				{
					uint64_t const ulow = t * supacksize;
					uint64_t const uhigh = std::min(ulow+supacksize,usize.back());
					assert ( uhigh != ulow );

					PP[uhigh] += PP[ulow];
				}

				std::vector < GenericWorkPackage < PackagePrefixSums > > VPPS(suthreads);
				for ( uint64_t t = 0; t < suthreads; ++t )
				{
					VPPS[t] = GenericWorkPackage < PackagePrefixSums >(0,PrefixSumsDispatcher_id,PackagePrefixSums(t,supacksize,&usize,&GRKM,&PP,&finsem));
					STP.enque(&VPPS[t]);
				}
				semWait(finsem,suthreads);

				std::vector<uint64_t> uranges;
				std::vector<uint64_t> urefk;
				uint64_t const uspan = (PP[usize.back()] + numthreads - 1)/numthreads;

				for ( uint64_t i = 0; i < numthreads; ++i )
				{
					uint64_t const utarget = i * uspan;
					uint64_t const * pp = std::lower_bound(PP.begin(),PP.begin()+usize.back(),utarget);
					uranges.push_back(pp-PP.begin());
					urefk.push_back(*pp);
				}
				uranges.push_back(usize.back());
				urefk.push_back(PP[usize.back()]);
				assert ( uranges[0] == 0 );

				#if 0
				for ( uint64_t i = 0; i < uranges.size(); ++i )
				{
					std::cerr << "uranges[.]=" << uranges[i]
						<< ((i > 0) ? (uranges[i]-uranges[i-1]) : 0)
						<< std::endl;
				}
				#endif

				#if 0
				{
					uint64_t s = 0;
					for ( uint64_t i = 0; i < usize.back(); ++i )
					{
						assert ( PP[i] == s );
						s += GRKM[i].high-GRKM[i].low;
					}
					assert ( PP[usize.back()] == s );
				}
				#endif

				PP.release();
				if ( verbose > 1 )
					std::cerr << "[V] computed range prefix sums in time " << lclock.formatTime(lclock.getElapsedSeconds()) << std::endl;

				// lookup positions for ranks on suffix array
				lclock.start();
				GRQKM.ensureSize(gc);
				std::vector<uint64_t> NUMREFKVALID(uranges.size());
				uint64_t const numlsa = uranges.size() ? (uranges.size()-1) : 0;
				std::vector < GenericWorkPackage < PackageLookupSA > > VLSA(numlsa);
				for ( uint64_t t = 1; t < uranges.size(); ++t )
				{
					VLSA[t-1] = GenericWorkPackage < PackageLookupSA >(0,LookupSADispatcher_id,PackageLookupSA(t,k,&GRKM,&GRQKM,&uranges,&urefk,&NUMREFKVALID,&index,&finsem));
					STP.enque(&VLSA[t-1]);
				}
				semWait(finsem,numlsa);

				// concatenate
				libmaus2::util::PrefixSums::prefixSums(NUMREFKVALID.begin(),NUMREFKVALID.end());
				TGRQKM.ensureSize(NUMREFKVALID.back()+2);
				uint64_t const numcomp = uranges.size() ? (uranges.size()-1) : 0;
				std::vector < GenericWorkPackage < PackageCompactPos > > VCO(numcomp);
				for ( uint64_t t = 1; t < uranges.size(); ++t )
				{
					VCO[t-1] = GenericWorkPackage < PackageCompactPos >(0,CompactPosDispatcher_id,PackageCompactPos(t,&GRQKM,&TGRQKM,&NUMREFKVALID,&urefk,&finsem));
					STP.enque(&VCO[t-1]);
				}
				semWait(finsem,numcomp);
				GRQKM.swap(TGRQKM);

				assert ( NUMREFKVALID.size() );
				assert ( GRQKM.size() >= NUMREFKVALID.back() );
				assert ( TGRQKM.size() >= NUMREFKVALID.back() );

				// sort
				libmaus2::sorting::InterleavedRadixSort::byteradixsortKeyBytes(
					GRQKM.begin(),GRQKM.begin()+NUMREFKVALID.back(),
					TGRQKM.begin(),TGRQKM.begin()+NUMREFKVALID.back(),
					numthreads,
					qfullkeybytes.begin(),qfullkeybytes.size(),STP
				);
				if ( qfullkeybytes.size() % 2 != 0 )
					GRQKM.swap(TGRQKM);
				#if 0
				std::cerr << "CHECKING" << std::endl;
				for ( uint64_t i = 1; i < NUMREFKVALID.back(); ++i )
					assert ( GRQKM[i-1] < GRQKM[i] );
				#endif

				GRQKM.ensureSize(NUMREFKVALID.back() + 2);

				GRQKM[NUMREFKVALID.back()] = QIff;
				GRQKM[NUMREFKVALID.back()+1] = QI0;
				if ( verbose > 1 )
					std::cerr << "[V] looked up and sorted positions in time " << lclock.formatTime(lclock.getElapsedSeconds()) << std::endl;

				lclock.start();
				libmaus2::util::PrefixSums::prefixSums(NUMQK.begin(),NUMQK.begin()+numpacks+1);
				GQKM.ensureSize(NUMQK[numpacks]+2);

				std::vector < GenericWorkPackage < PackageConcatQueryKmers > > VCQK(numpacks);
				for ( uint64_t t = 0; t < numpacks; ++t )
				{
					VCQK[t] = GenericWorkPackage < PackageConcatQueryKmers >(0,ConcatQueryKmersDispatcher_id,PackageConcatQueryKmers(t,&QKM,&GQKM,&NUMQK,&finsem));
					STP.enque(&VCQK[t]);
				}
				semWait(finsem,numpacks);
				assert ( ! finsem.trywait() );
				if ( verbose > 1 )
					std::cerr << "[V] concatenated query kmers in time " << lclock.formatTime(lclock.getElapsedSeconds()) << std::endl;

				lclock.start();
				TGQKM.ensureSize(NUMQK[numpacks]+2);
				// we sort the query kmers by the code field only, because readid and rpos are produced already in order and radix sort is stable
				libmaus2::sorting::InterleavedRadixSort::byteradixsortKeyBytes(
					GQKM.begin(),GQKM.begin()+NUMQK[numpacks],
					TGQKM.begin(),TGQKM.begin()+NUMQK[numpacks],
					numthreads,
					keybytes.begin(),numkeybytes,STP);
				if ( verbose > 1 )
					std::cerr << "[V] sorted query kmers in time " << lclock.formatTime(lclock.getElapsedSeconds()) << std::endl;
				if ( numkeybytes % 2 == 1 )
					std::swap(GQKM,TGQKM);

				#if 0
				std::cerr << "checking q seeds" << std::endl;
				for ( uint64_t i = 1; i < NUMQK[numpacks]; ++i )
				{
					assert (
						GQKM[i-1].code < GQKM[i].code
						||
						(
							GQKM[i-1].code == GQKM[i].code &&
							GQKM[i-1].readid < GQKM[i].readid
						)
						||
						(
						GQKM[i-1].code == GQKM[i].code &&
						GQKM[i-1].readid == GQKM[i].readid &&
						GQKM[i-1].rpos < GQKM[i].rpos
						)
					);
				}
				#endif

				GQKM[NUMQK[numpacks]] = QIff;
				GQKM[NUMQK[numpacks]+1] = QI0;

				#if 0
				for ( uint64_t i = 1; i < NUMQK[numpacks]; ++i )
					assert ( GQKM[i-1].code <= GQKM[i].code );
				for ( uint64_t i = 1; i < NUMREFKVALID.back(); ++i )
					assert ( GRQKM[i-1].code <= GRQKM[i].code );
				#endif

				// decode ref sequences for this index
				lclock.start();
				std::vector<libmaus2::fastx::FastaBPDecoderIdentity::SequenceMeta> Vseqmeta;
				libmaus2::autoarray::AutoArray<char,libmaus2::autoarray::alloc_type_c> Aseq;
				index.loadBP(STP,Vseqmeta,Aseq,numthreads);
				index.fillDazzlerDB(refdb,AREFHITREADS,Vseqmeta,Aseq.begin());
				if ( verbose > 1 )
					std::cerr << "[V] decode ref seqs in time " << lclock.formatTime(lclock.getElapsedSeconds()) << std::endl;

				#if 0
				std::cerr << "checking seeds" << std::endl;
				// check seed correctness by comparing reference and query text
				QueryKMer * ref_a = GRQKM.begin();
				QueryKMer * ref_e = ref_a + NUMREFKVALID.back();
				QueryKMer * q_a = GQKM.begin();
				QueryKMer * q_e = GQKM.begin() + NUMQK[numpacks];

				while ( ref_a != ref_e && q_a != q_e )
				{
					if ( ref_a->code < q_a->code )
						++ref_a;
					else if ( q_a->code < ref_a->code )
						++q_a;
					else
					{
						uint64_t const code = q_a->code;
						QueryKMer * ref_s = ref_a;
						QueryKMer * q_s = q_a;

						while ( ref_a != ref_e && ref_a->code == code )
							++ref_a;
						while ( q_a != q_e && q_a->code == code )
							++q_a;

						for ( QueryKMer * ref_c = ref_s; ref_c != ref_a; ++ref_c )
							for ( QueryKMer * q_c = q_s; q_c != q_a; ++ q_c )
							{
								assert ( ref_c->code == q_c->code );

								char const * aseq = reinterpret_cast<char const *>(readsdb.bases) + readsdb.reads[q_c->readid].boff;
								char const * bseq = reinterpret_cast<char const *>(refdb.bases)   + refdb.reads[ref_c->readid].boff;

								assert ( aseq[-1] == 4 );
								assert ( aseq[readsdb.reads[q_c->readid].rlen] == 4 );
								assert ( bseq[-1] == 4 );
								assert ( bseq[refdb.reads[ref_c->readid].rlen] == 4 );

								// std::cerr << std::string(80,'-') << std::endl;

								for ( uint64_t z = 0; z < k; ++z )
								{
									// std::cerr << (int)aseq [ SP.apos + z ] << "\t" << (int)bseq [ SP.apos - SP.diag + z ] << std::endl;

									assert (
										aseq [ q_c->rpos + z ] ==
										bseq [ ref_c->rpos + z ]
									);
								}

							}
					}
				}
				#endif
				GRKM.release();
				TGRKM.release();
				TGQKM.release();
				TGRQKM.release();

				if ( NUMQK[numpacks] && NUMREFKVALID.back() )
				{
					lclock.start();
					QueryKMer * Qindex = GRQKM.take();
					Match_Filter(&readsdb, &refdb, GQKM.begin(), NUMQK[numpacks], Qindex, NUMREFKVALID.back(), (index_i==1) /* comp */, (index_i==0) /* start */);
					indexyield += 1;
					if ( verbose > 1 )
						std::cerr << "[V] ran Match_Filter in time " << lclock.formatTime(lclock.getElapsedSeconds()) << std::endl;
				}
				else
				{
					GRQKM.release();
				}

				GQKM.release();
			}

			lclock.start();
			std::vector<libmaus2::fastx::FastaBPDecoderIdentity::SequenceMeta> Vseqmeta;
			libmaus2::autoarray::AutoArray<char,libmaus2::autoarray::alloc_type_c> Aseq;
			forwindex.loadBP(STP,Vseqmeta,Aseq,numthreads);
			forwindex.fillDazzlerDB(refdb,AREFHITREADS,Vseqmeta,Aseq.begin());
			if ( verbose > 1 )
				std::cerr << "[V] decode ref seqs in time " << lclock.formatTime(lclock.getElapsedSeconds()) << std::endl;

			std::string const tmpprefix = libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname);
			libmaus2::util::WriteableString Wprefix(tmpprefix + "_reads");
			libmaus2::util::WriteableString Wref("ref");

			std::vector<std::string> Voutfn;

			if ( indexyield )
			{
				Reporter(Wprefix.A.begin(), &readsdb, Wref.A.begin(), &refdb, aspec, FLAG_DOB /* mflag */);

				// for ( uint64_t i = 1; i <= NTHREADS; ++i )
				for ( uint64_t i = 1; i <= numthreads; ++i )
				{
					std::ostringstream fnostr;
					#if 0
					fnostr
						<< Wref.A.begin()
						<< "."
						<< Wprefix.A.begin()
						<< ".R" << i << ".las";
					#else
					fnostr
						<< "/tmp/"
						<< Wref.A.begin()
						<< "."
						<< Wprefix.A.begin()
						<< ".R" << i << ".las";
					#endif
					Voutfn.push_back(fnostr.str());
				}
			}

			libmaus2::dazzler::align::SimpleOverlapVectorParser::unique_ptr_type PSOP(
				new libmaus2::dazzler::align::SimpleOverlapVectorParser(
					Voutfn,1024*1024 /* bufsize */,
					libmaus2::dazzler::align::OverlapParser::overlapparser_do_not_split_b
				)
			);

			int64_t prevmark = -1;
			int64_t maxbread = -1;

			while ( PSOP->parseNextBlock() )
			{
				libmaus2::dazzler::align::OverlapData & OVLdata = PSOP->getData();

				uint64_t rlow = 0;
				uint64_t rsum = 0;
				uint64_t linesperthread = (OVLdata.size() + numthreads - 1) / numthreads;
				std::vector<uint64_t> ointv;
				ointv.push_back(0);

				while ( rlow < OVLdata.size() )
				{
					// reference b read id
					int64_t const refb = libmaus2::dazzler::align::OverlapData::getBRead(OVLdata.getData(rlow).first);
					uint64_t rhigh = rlow+1;

					// look for end of b read range
					while ( rhigh < OVLdata.size() && libmaus2::dazzler::align::OverlapData::getBRead(OVLdata.getData(rhigh).first) == refb )
						++rhigh;

					// update sum
					rsum += (rhigh-rlow);

					// if we have enough
					if ( rsum >= linesperthread )
					{
						rsum = 0;
						ointv.push_back(rhigh);
					}

					rlow = rhigh;
				}

				if ( rsum )
				{
					rsum = 0;
					ointv.push_back(rlow);
				}

				for ( uint64_t i = 1; i < ointv.size(); ++i )
					assert ( ointv[i] > ointv[i-1] );

				assert ( ointv.size() > 0 );
				assert ( ointv.size()-1 <= numthreads );
				assert ( ointv.back() == OVLdata.size() );

				uint64_t const numointv = ointv.size() ? (ointv.size()-1) : 0;
				std::vector < GenericWorkPackage < PackageEncodeBam > > VEB(numointv);
				for ( uint64_t zr = 1; zr < ointv.size(); ++zr )
				{
					VEB[zr-1] = GenericWorkPackage < PackageEncodeBam >(
						0,EncodeBamDispatcher_id,
						PackageEncodeBam(zr,&ointv,&Acontexts,&OVLdata,prevmark,&Areadnames,&Pbamheader,&O,&refdb,&readsdb,&Abam,&Obam,oobam,&finsem)
					);
					STP.enque(&VEB[zr-1]);
				}
				semWait(finsem,numointv);

				std::vector<std::pair<uint8_t *,uint8_t *> > Vlinfrag;
				for ( uint64_t zr = 1; zr < ointv.size(); ++zr )
				{
					uint64_t const tid = zr-1;
					LASToBamContext & context = *(Acontexts[tid]);
					libmaus2::bambam::parallel::FragmentAlignmentBufferFragment & fragment = context.fragment;
					fragment.getLinearOutputFragments(libmaus2::lz::BgzfConstants::getBgzfMaxBlockSize(),Vlinfrag);
				}


				std::vector < libmaus2::lz::BgzfDeflateZStreamBaseFlushInfo > VBinfo(Vlinfrag.size());
				std::vector < GenericWorkPackage < PackageCompressBgzf > > DVCB(Vlinfrag.size());
				for ( uint64_t i = 0; i < Vlinfrag.size(); ++i )
				{
					DVCB[i] = GenericWorkPackage < PackageCompressBgzf >(0/*prio*/,CompressBgzfDispatcher_id,PackageCompressBgzf(&zstreambasefreelist,Vlinfrag[i],bgzffreelist.get(),&VBinfo[i],&finsem));
					STP.enque(&DVCB[i]);
				}
				semWait(finsem,DVCB.size());

				for ( uint64_t i = 0; i < DVCB.size(); ++i )
				{
					uint8_t const * u = DVCB[i].package.B->outbuf.begin();
					std::cout.write(reinterpret_cast<char const *>(u),DVCB[i].package.info->getCompressedSize());
					bgzffreelist.put(DVCB[i].package.B);
				}

				std::cerr << "[V] cap " << zstreambasefreelist.capacity() << " " << bgzffreelist.capacity() << std::endl;

				if ( OVLdata.size() )
				{
					prevmark = libmaus2::dazzler::align::OverlapData::getBRead(OVLdata.getData(OVLdata.size()-1).first);
					maxbread = std::max(maxbread,prevmark);
				}

				if ( verbose > 1 )
					std::cerr << "[V] processed " << OVLdata.size() << " alignments " << libmaus2::util::MemUsage() << std::endl;
			}

			std::cerr << "[V] maxbread=" << maxbread << " numreads=" << numreads << std::endl;

			for ( uint64_t i = maxbread + 1; i < numreads; ++i )
			{
				std::cerr << "[V] read " << i << " missing in output" << std::endl;
			}

			running = havedata;

			totalreads += numreads;

			if ( verbose > 0 )
			{
				double const rcl = roundclock.getElapsedSeconds();
				double const acccl = accclock.getElapsedSeconds();
				double const roundalpersec = numreads / rcl;
				double const accalpersec = totalreads / acccl;

				std::cerr << "[V]"
					<< " processed " << numreads << " reads in time " << roundclock.formatTime(rcl)
					<< " " << roundalpersec << " al/s"
					<< " total " << totalreads << " time " << accclock.formatTime(acccl)
					<< " " << accalpersec << " al/s"
					<< " " << accalpersec*60*60 << " al/h"
					<< std::endl;
			}

			for ( uint64_t i = 0; i < Voutfn.size(); ++i )
				libmaus2::aio::FileRemoval::removeFile(Voutfn[i]);
		}

		// bgzfout->flush();
		bgzfout->addEOFBlock();
		bgzfout.reset();

		Free_Align_Spec(aspec);

		STP.terminate();

		return EXIT_SUCCESS;
	}
	catch(...)
	{
		STP.terminate();
		throw;
	}
}

#include <config.h>

std::string getUsage(libmaus2::util::ArgParser const & arg)
{
	std::ostringstream ostr;

	ostr << "usage: " << arg.progname << " [<parameters>] ref.fasta <reads.fasta" << std::endl;
	ostr << "\n";
	ostr << "parameters:\n";
	ostr << " -k: seed length (default 20)\n";
	ostr << " -K: kmer cache word length (default 12)\n";
	#if defined(_OPENMP)
	ostr << " -p: number of threads (defaults to number of cores on machine)\n";
	#else
	ostr << " -p: no effect (compiled without support for OpenMP parallelism)\n";
	#endif
	ostr << " -v: verbosity level (default: 1)\n";
	ostr << " -i: maximum number of read bases per input block (default 256m)\n";
	ostr << " -M: damapper memory limit\n";
	ostr << " --sasamplingrate: SA sampling rate (default 4)\n";
	ostr << " --bwtconstrmem: memory used to construct BWT (default 3/4 of machine's memory)\n";
	ostr << " -T: prefix for temporary files used during index construction\n";
	ostr << " -Q: file name of index\n";
	ostr << " -S: storage strategy for non primary alignments (none, soft, hard)\n";
	ostr << " -z: output BAM compression level (zlib default)\n";
	ostr << " -I: input format (fasta (default) or bam)\n";

	return ostr.str();
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgParser const arg(argc,argv);

		if ( arg.argPresent("h") || arg.argPresent("help") )
		{
			std::cerr << getUsage(arg);
			return EXIT_SUCCESS;
		}
		else if ( arg.argPresent("version") )
		{
			std::cerr << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << std::endl;
			return EXIT_SUCCESS;
		}
		else if ( arg.size() < 1 )
		{
			std::cerr << getUsage(arg);
			return EXIT_FAILURE;
		}

		return damapper_bwt(arg);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
