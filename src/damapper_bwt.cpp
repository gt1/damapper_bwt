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

#include <libmaus2/aio/SynchronousGenericInput.hpp>
#include <libmaus2/bambam/BamAlignmentEncoderBase.hpp>
#include <libmaus2/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus2/bambam/BamHeader.hpp>
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

	libmaus2::bambam::BamHeader::unique_ptr_type getBamHeader()
	{
		std::ostringstream samheaderstr;
		samheaderstr << "@HD\tVN:1.5\tSO:unknown\n";
		for ( uint64_t i = 0; i < PFAI->size(); ++i )
		{
			libmaus2::fastx::FastAIndexEntry const & entry = (*PFAI)[i];
			samheaderstr << "@SQ\tSN:" << entry.name << "\tLN:" << entry.length << "\n";
		}
		// std::cerr << samheaderstr.str();
		libmaus2::bambam::BamHeader::unique_ptr_type Pbamheader(new libmaus2::bambam::BamHeader(samheaderstr.str()));

		return UNIQUE_PTR_MOVE(Pbamheader);
	}

	void loadBP(std::vector<libmaus2::fastx::FastaBPDecoderIdentity::SequenceMeta> & Vseqmeta, libmaus2::autoarray::AutoArray<char,libmaus2::autoarray::alloc_type_c> & Aseq, uint64_t const numthreads)
	{
		libmaus2::aio::InputStreamInstance::unique_ptr_type membpISI(new libmaus2::aio::InputStreamInstance(membpname));
		libmaus2::fastx::FastaBPDecoderIdentity fabpdec(*membpISI);
		fabpdec.decodeSequencesParallel(membpname,numthreads,Aseq,Vseqmeta,true /* map */,4 /* pad symbol */,false /* addrc */);
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
		refdb.all = 0;
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

int damapper_bwt(libmaus2::util::ArgParser const & arg)
{

	// default k
	unsigned int const defk = 20;
	unsigned int const k = arg.argPresent("k") ? arg.getUnsignedNumericArg<uint64_t>("k") : defk;

	// verbosity
	unsigned int const defv = 1;
	unsigned int const verbose = arg.argPresent("v") ? arg.getUnsignedNumericArg<uint64_t>("v") : defv;

	// maximum number of input bases per block
	uint64_t const defo = 64*1024*1024;
	uint64_t const maxo = arg.argPresent("i") ? std::max(static_cast<uint64_t>(1),arg.getUnsignedNumericArg<uint64_t>("i")) : defo;

	// length for K-mer cache
	unsigned int const defcachek = 12;
	unsigned int const cache_k = arg.uniqueArgPresent("K") ? arg.getUnsignedNumericArg<uint64_t>("K") : defcachek;

	MEM_PHYSICAL = libmaus2::util::MemoryStatistics::getPhysicalMemory();
	MEM_LIMIT = MEM_PHYSICAL;

	uint64_t const defconstrsasamplingrate = 32;
	uint64_t const constrsasamplingrate = arg.argPresent("sasamplingrate") ? std::max(static_cast<uint64_t>(1),arg.getUnsignedNumericArg<uint64_t>("sasamplingrate")) : defconstrsasamplingrate;
	uint64_t const deficonstrsasamplingrate = 32;
	uint64_t const iconstrsasamplingrate = arg.argPresent("isasamplingrate") ? std::max(static_cast<uint64_t>(1),arg.getUnsignedNumericArg<uint64_t>("isasamplingrate")) : deficonstrsasamplingrate;
	uint64_t const defaultbwtconstrmem = (MEM_PHYSICAL*3)/4;
	uint64_t const bwtconstrmem = arg.argPresent("bwtconstrmem") ? std::max(static_cast<uint64_t>(1),arg.getUnsignedNumericArg<uint64_t>("bwtconstrmem")) : defaultbwtconstrmem;

	#if defined(_OPENMP)
	uint64_t const defnumproc = libmaus2::parallel::NumCpus::getNumLogicalProcessors();
	uint64_t const numthreads = arg.argPresent("p") ? std::max(static_cast<uint64_t>(1),arg.getUnsignedNumericArg<uint64_t>("p")) : defnumproc;
	#else
	uint64_t const numthreads = 1;
	#endif

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
	combISI.reset();

	libmaus2::bambam::BamHeader::unique_ptr_type Pbamheader(forwindex.getBamHeader());

	std::vector < DNAIndex * > indexes;
	indexes.push_back(&forwindex);
	indexes.push_back(&rcindex);

	float freq[4] = { 0.25, 0.25, 0.25, 0.25 };
	int64_t const tspace = 100;
	Align_Spec * aspec = New_Align_Spec(0.85, tspace, &freq[0]);
	Set_Filter_Params(k, 0, numthreads);

	libmaus2::fastx::LineBufferFastAReader LBFA(std::cin);
	bool running = true;

	libmaus2::timing::RealTimeClock roundclock;
	libmaus2::timing::RealTimeClock lclock;

	// read names
	libmaus2::autoarray::AutoArray<char> Areadnames;
	// read base data
	libmaus2::autoarray::AutoArray<char> Areaddata;
	// read meta data
	libmaus2::autoarray::AutoArray< libmaus2::fastx::LineBufferFastAReader::ReadMeta > O;

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

	libmaus2::lz::BgzfDeflate<std::ostream>::unique_ptr_type bgzfout(new libmaus2::lz::BgzfDeflate<std::ostream>(std::cout,Z_NO_COMPRESSION));
	Pbamheader->serialise(*bgzfout);

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

	libmaus2::autoarray::AutoArray<LASToBamContext::unique_ptr_type> Acontexts(numthreads);
	libmaus2::dazzler::align::RefMapEntryVector refmap;
	for ( uint64_t i = 0; i < Pbamheader->getNumRef(); ++i )
		refmap.push_back(libmaus2::dazzler::align::RefMapEntry(i,0));
	for ( uint64_t i = 0; i < numthreads; ++i )
	{
		LASToBamContext::unique_ptr_type Tptr(new LASToBamContext(tspace,true /* mdnm */,libmaus2::dazzler::align::LASToBamConverterBase::supplementary_seq_strategy_none,std::string(),refmap));
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

		lclock.start();
		// load read block
		while ( o_data < maxo && (havedata = LBFA.getNext(Areadnames,o_name,Areaddata,o_data,rlen,termdata,mapdata,paddata,padsym)) )
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
		readsdb.all = 0;
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
			#if defined(_OPENMP)
			#pragma omp parallel for num_threads(numthreads) schedule(dynamic,1)
			#endif
			for ( uint64_t t = 0; t < numpacks; ++t )
			{
				uint64_t const low = t * readsperpack;
				uint64_t const high = std::min(low+readsperpack,numreads);
				libmaus2::autoarray::AutoArray < RefKMer > & LRKM = *RKM[t];
				libmaus2::autoarray::AutoArray < QueryKMer > & LQKM = *QKM[t];

				uint64_t rkmo = 0;
				uint64_t qkmo = 0;

				uint64_t const kmask = libmaus2::math::lowbits(2*k);

				for ( uint64_t zz = low; zz < high; ++zz )
				{
					//char const * const cmapped = Areaddata.begin() + O[zz].dataoff;
					char * const cmapped = Areaddata.begin() + O[zz].dataoff;
					uint64_t const m = O[zz].len;

					uint64_t const numk = (m >= k) ? m-k+1 : 0;

					if ( numk && k )
					{
						uint64_t v = 0;
						uint64_t e = 0;
						char const * c = cmapped;

						for ( uint64_t i = 0; i < k-1; ++i )
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
								std::pair<uint64_t,uint64_t> const P = index.PKcache->search(c-k,k);

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

				NUMRK[t] = rkmo;
				NUMQK[t] = qkmo;
			}
			if ( verbose > 1 )
				std::cerr << "[V] computed SA ranges in time " << lclock.formatTime(lclock.getElapsedSeconds()) << std::endl;

			// concatenate reference kmer sequences
			lclock.start();
			libmaus2::util::PrefixSums::prefixSums(NUMRK.begin(),NUMRK.begin()+numpacks+1);
			GRKM.ensureSize(NUMRK[numpacks]);
			#if defined(_OPENMP)
			#pragma omp parallel for num_threads(numthreads) schedule(dynamic,1)
			#endif
			for ( uint64_t t = 0; t < numpacks; ++t )
			{
				RefKMer * O = GRKM.begin() + NUMRK[t];
				RefKMer * Oe = GRKM.begin() + NUMRK[t+1];
				RefKMer * I = RKM[t]->begin();

				while ( O != Oe )
					*(O++) = *(I++);

				RKM[t]->resize(0);
			}
			if ( verbose > 1 )
				std::cerr << "[V] concatenated SA ranges in time " << lclock.formatTime(lclock.getElapsedSeconds()) << std::endl;

			// sort suffix ranges by code field (so we can unify the ranges)
			lclock.start();
			TGRKM.ensureSize(NUMRK[numpacks]);
			libmaus2::sorting::InterleavedRadixSort::byteradixsortKeyBytes(
				GRKM.begin(),GRKM.begin()+NUMRK[numpacks],
				TGRKM.begin(),TGRKM.begin()+NUMRK[numpacks],
				numthreads,
				keybytes.begin(),numkeybytes);
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

			#if defined(_OPENMP)
			#pragma omp parallel for num_threads(numthreads)
			#endif
			for ( uint64_t t = 0; t < numthreads; ++t )
			{
				std::fill(
					Ahist[t]->begin(),
					Ahist[t]->end(),
					0ull
				);
			}

			// unified package sizes
			std::vector<uint64_t> usize(ubounds.size());
			#if defined(_OPENMP)
			#pragma omp parallel for num_threads(numthreads)
			#endif
			for ( uint64_t t = 1; t < ubounds.size(); ++t )
			{
				libmaus2::autoarray::AutoArray<uint64_t> & thist = *(Ahist[t-1].get());

				// lower bound
				uint64_t ulow = ubounds[t-1];
				// upper bound
				uint64_t const utop = ubounds[t];

				// output pointer
				uint64_t uout = ulow;
				while ( ulow < utop )
				{
					// go to end of code
					uint64_t uhigh = ulow+1;
					while ( uhigh < utop && GRKM[uhigh].code == GRKM[ulow].code )
						++uhigh;

					// single output
					GRKM[uout++] = GRKM[ulow];

					uint64_t const s = GRKM[ulow].high-GRKM[ulow].low;
					if ( s < histthres )
						thist [ s ] += 1;
					else
						thist [ histthres ] += s;

					ulow = uhigh;
				}

				// size of package
				usize[t-1] = uout - ubounds[t-1];
			}
			// compute prefix sums
			libmaus2::util::PrefixSums::prefixSums(usize.begin(),usize.end());

			// merge histograms
			uint64_t const tperthread = (histmult + numthreads - 1)/numthreads;
			#if defined(_OPENMP)
			#pragma omp parallel for num_threads(numthreads)
			#endif
			for ( uint64_t t = 0; t < numthreads; ++t )
			{
				uint64_t const tlow = t * tperthread;
				uint64_t const thigh = std::min ( tlow + tperthread, histmult );

				uint64_t * const T = Ahist[0]->begin();

				for ( uint64_t j = 1; j < numthreads; ++j )
				{
					uint64_t * const S = Ahist[j]->begin();

					for ( uint64_t i = tlow; i < thigh; ++i )
						T[i] += S[i];
				}
			}

			#if 0
			for ( uint64_t i = 0; i < histmult; ++i )
				if ( (*(Ahist[0]))[i] )
				{
					std::cerr << "[H] " << i << " " << (*(Ahist[0]))[i] << std::endl;
				}
			#endif

			// concatenate unique
			uint64_t gc = 0;
			libmaus2::parallel::PosixSpinLock gclock;
			#if defined(_OPENMP)
			#pragma omp parallel for num_threads(numthreads)
			#endif
			for ( uint64_t t = 1; t < usize.size(); ++t )
			{
				RefKMer * I = GRKM.begin() + ubounds[t-1];
				RefKMer * O = TGRKM.begin() + usize[t-1];
				RefKMer * Oe = TGRKM.begin() + usize[t];

				uint64_t c = 0;
				while ( O != Oe )
				{
					c += I->high-I->low;
					*(O++) = *(I++);
				}
				gclock.lock();
				gc += c;
				gclock.unlock();

				#if 0
				uint64_t cc = 0;
				for ( RefKMer * C = TGRKM.begin() + usize[t-1]; C != TGRKM.begin() + usize[t]; ++C )
					cc += C->high-C->low;
				assert ( cc == c );
				#endif
			}
			GRKM.swap(TGRKM);
			if ( verbose > 1 )
				std::cerr << "[V] uniquified in time " << lclock.formatTime(lclock.getElapsedSeconds()) << " num out " << usize.back() << " sum " << gc << std::endl;

			lclock.start();
			libmaus2::autoarray::AutoArray<uint64_t,libmaus2::autoarray::alloc_type_c> PP(usize.back()+1,false);
			uint64_t const supacksize = (usize.back() + numthreads - 1)/numthreads;
			uint64_t const suthreads = supacksize ? ((usize.back() + supacksize - 1) /supacksize) : 0;
			PP[0] = 0;

			#if defined(_OPENMP)
			#pragma omp parallel for num_threads(suthreads)
			#endif
			for ( uint64_t t = 0; t < suthreads; ++t )
			{
				uint64_t const ulow = t * supacksize;
				uint64_t const uhigh = std::min(ulow+supacksize,usize.back());
				assert ( uhigh != ulow );

				RefKMer * I  = GRKM.begin() + ulow;
				RefKMer * Ie = GRKM.begin() + uhigh;

				uint64_t s = 0;
				for ( ; I != Ie ; ++I )
					s += I->high-I->low;

				PP[uhigh] = s;
			}

			for ( uint64_t t = 0; t < suthreads; ++t )
			{
				uint64_t const ulow = t * supacksize;
				uint64_t const uhigh = std::min(ulow+supacksize,usize.back());
				assert ( uhigh != ulow );

				PP[uhigh] += PP[ulow];
			}

			#if defined(_OPENMP)
			#pragma omp parallel for num_threads(suthreads)
			#endif
			for ( uint64_t t = 0; t < suthreads; ++t )
			{
				uint64_t const ulow = t * supacksize;
				uint64_t const uhigh = std::min(ulow+supacksize,usize.back());
				assert ( uhigh != ulow );

				RefKMer * I  = GRKM.begin() + ulow;
				RefKMer * Ie = GRKM.begin() + uhigh;

				uint64_t * p = PP.begin() + ulow;

				uint64_t s = *p;
				for ( ; I != Ie ; ++I )
				{
					*(p++) = s;
					s += (I->high-I->low);
				}
			}

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
			#if defined(_OPENMP)
			#pragma omp parallel for num_threads(numthreads)
			#endif
			for ( uint64_t t = 1; t < uranges.size(); ++t )
			{
				RefKMer * I  = GRKM.begin() + uranges[t-1];
				RefKMer * Ie = GRKM.begin() + uranges[t];
				QueryKMer * O = GRQKM.begin() + urefk[t-1];
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
						uint64_t const p = index.lookupSA(jr);
						std::pair<uint64_t,uint64_t> coord = (*(index.PCC))[p]; // Pindex->mapCoordinates((*SSA)[jr]);
						if ( index.Pindex->valid(coord,k) )
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
				assert ( c == urefk[t]-urefk[t-1] );
				// assert ( O == GRQKM.begin() + NUMREFK[t] );
				NUMREFKVALID[t-1] = numvalid;
			}

			// concatenate
			libmaus2::util::PrefixSums::prefixSums(NUMREFKVALID.begin(),NUMREFKVALID.end());
			TGRQKM.ensureSize(NUMREFKVALID.back()+2);
			#if defined(_OPENMP)
			#pragma omp parallel for num_threads(numthreads)
			#endif
			for ( uint64_t t = 1; t < uranges.size(); ++t )
			{
				QueryKMer * I = GRQKM.begin() + urefk[t-1];
				QueryKMer * O = TGRQKM.begin() + NUMREFKVALID[t-1];
				QueryKMer * Oe = TGRQKM.begin() + NUMREFKVALID[t];

				while ( O != Oe )
					*(O++) = *(I++);
			}
			GRQKM.swap(TGRQKM);

			assert ( NUMREFKVALID.size() );
			assert ( GRQKM.size() >= NUMREFKVALID.back() );
			assert ( TGRQKM.size() >= NUMREFKVALID.back() );

			// sort
			libmaus2::sorting::InterleavedRadixSort::byteradixsortKeyBytes(
				GRQKM.begin(),GRQKM.begin()+NUMREFKVALID.back(),
				TGRQKM.begin(),TGRQKM.begin()+NUMREFKVALID.back(),
				numthreads,
				qfullkeybytes.begin(),qfullkeybytes.size()
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

			#if defined(_OPENMP)
			#pragma omp parallel for num_threads(numthreads) schedule(dynamic,1)
			#endif
			for ( uint64_t t = 0; t < numpacks; ++t )
			{
				libmaus2::autoarray::AutoArray < QueryKMer > & LQKM = *QKM[t];

				QueryKMer const * I = LQKM.begin();
				QueryKMer * O = GQKM.begin() + NUMQK[t];
				QueryKMer * Oe = GQKM.begin() + NUMQK[t+1];

				while ( O != Oe )
					*(O++) = *(I++);

				QKM[t]->resize(0);
			}
			if ( verbose > 1 )
				std::cerr << "[V] concatenated query kmers in time " << lclock.formatTime(lclock.getElapsedSeconds()) << std::endl;

			lclock.start();
			TGQKM.ensureSize(NUMQK[numpacks]+2);
			// we sort the query kmers by the code field only, because readid and rpos are produced already in order and radix sort is stable
			libmaus2::sorting::InterleavedRadixSort::byteradixsortKeyBytes(
				GQKM.begin(),GQKM.begin()+NUMQK[numpacks],
				TGQKM.begin(),TGQKM.begin()+NUMQK[numpacks],
				numthreads,
				keybytes.begin(),numkeybytes);
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
			index.loadBP(Vseqmeta,Aseq,numthreads);
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
				QueryKMer * Qindex = GRQKM.take();
				Match_Filter(&readsdb, &refdb, GQKM.begin(), NUMQK[numpacks], Qindex, NUMREFKVALID.back(), (index_i==1) /* comp */, (index_i==0) /* start */);
				indexyield += 1;
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
		forwindex.loadBP(Vseqmeta,Aseq,numthreads);
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
				fnostr
					<< Wref.A.begin()
					<< "."
					<< Wprefix.A.begin()
					<< ".C" << i << ".las";
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

			#if defined(_OPENMP)
			#pragma omp parallel for num_threads(numthreads)
			#endif
			for ( uint64_t zr = 1; zr < ointv.size(); ++zr )
			{
				uint64_t tid = zr-1;
				uint64_t const rlow = ointv[zr-1];
				uint64_t const rhigh = ointv[zr];
				LASToBamContext & context = *(Acontexts[tid]);
				context.fragment.reset();

				// last read id for previous block
				int64_t const prevb = rlow ? static_cast<int64_t>(libmaus2::dazzler::align::OverlapData::getBRead(OVLdata.getData(rlow-1).first)) : prevmark;
				// next expected read
				int64_t nextexpt = prevb+1;

				//std::cerr << "tid " << tid << " " << nextexpt << std::endl;

				// process in batches of equal b read id
				uint64_t ilow = rlow;
				while ( ilow < rhigh )
				{
					uint64_t ihigh = ilow+1;
					int64_t const refbread = libmaus2::dazzler::align::OverlapData::getBRead(OVLdata.getData(ilow).first);
					while ( ihigh < rhigh && libmaus2::dazzler::align::OverlapData::getBRead(OVLdata.getData(ihigh).first) == refbread )
						++ihigh;

					while ( nextexpt < refbread )
					{
						char const * rname = Areadnames.begin() + O[nextexpt].nameoff;
						char const * rdata = reinterpret_cast<char const *>(readsdb.bases) + readsdb.reads[nextexpt].boff;
						uint64_t const readlen = readsdb.reads[nextexpt].rlen;
						assert ( rdata[-1] == 4 );
						assert ( rdata[readlen] == 4 );

						context.converter.convertUnmapped(rdata,readlen,rname,context.fragment,*Pbamheader);

						// std::cerr << "[D] read " << nextexpt << " " << rname << " was not aligned" << " ref " << refbread << std::endl;

						++nextexpt;
					}

					assert ( nextexpt == refbread );

					uint64_t numchains = 0;
					for ( uint64_t z = ilow; z < ihigh; ++z )
					{
						std::pair<uint8_t const *, uint8_t const *> const P = OVLdata.getData(z);
						bool const isStart = libmaus2::dazzler::align::OverlapData::getStartFlag(P.first);
						if ( isStart )
							numchains++;
					}

					uint64_t il = ilow;
					uint64_t chainid = 0;
					while ( il < ihigh )
					{
						uint64_t ih = il+1;
						while (
							ih < ihigh
							&&
							(!libmaus2::dazzler::align::OverlapData::getStartFlag(OVLdata.getData(ih).first))
						)
							++ih;

						uint64_t const lchainid = chainid++;
						bool const secondary = (lchainid > 0);

						for ( uint64_t z = il; z < ih; ++z )
						{
							std::pair<uint8_t const *, uint8_t const *> const P = OVLdata.getData(z);

							int64_t const aread = libmaus2::dazzler::align::OverlapData::getARead(P.first);
							int64_t const bread = libmaus2::dazzler::align::OverlapData::getBRead(P.first);
							int64_t const flags = libmaus2::dazzler::align::OverlapData::getFlags(P.first);
							bool const inverse = libmaus2::dazzler::align::OverlapData::getInverseFlag(P.first);

							bool const supplementary = (z != il);

							// bool const primary = libmaus2::dazzler::align::OverlapData::getPrimaryFlag(P.first);
							uint64_t const readlen = readsdb.reads[bread].rlen;

							if ( z==ilow )
							{
								// compute reverse complement
								context.ARC.ensureSize(readlen + 2);

								char * const ra = context.ARC.begin();
								char * rp = ra + readlen + 2;
								char const * src = reinterpret_cast<char const *>(readsdb.bases) + readsdb.reads[bread].boff - 1;

								*(--rp)	= (*(src++));
								while ( rp != ra+1 )
									*(--rp)	= (*(src++)) ^ 3;
								*(--rp)	= (*(src++));
							}

							libmaus2::dazzler::align::LASToBamConverterBase::AuxTagIntegerAddRequest req_i("ci",lchainid);
							libmaus2::dazzler::align::LASToBamConverterBase::AuxTagIntegerAddRequest req_n("cn",numchains);
							libmaus2::dazzler::align::LASToBamConverterBase::AuxTagIntegerAddRequest req_j("cj",z-il);
							libmaus2::dazzler::align::LASToBamConverterBase::AuxTagIntegerAddRequest req_l("cl",ih-il);

							libmaus2::dazzler::align::LASToBamConverterBase::AuxTagAddRequest const * reqs[] =
							{
								&req_i,
								&req_n,
								&req_j,
								&req_l
							};

							libmaus2::dazzler::align::LASToBamConverterBase::AuxTagAddRequest const ** aux_a = &reqs[0];
							libmaus2::dazzler::align::LASToBamConverterBase::AuxTagAddRequest const ** aux_e = &reqs[sizeof(reqs)/sizeof(reqs[0])];

							try
							{
								if ( ! inverse )
								{
									context.converter.convert(
										P.first,
										reinterpret_cast<char const *>(refdb.bases) + refdb.reads[aread].boff,
										refdb.reads[aread].rlen,
										reinterpret_cast<char const *>(readsdb.bases) + readsdb.reads[bread].boff,
										readsdb.reads[bread].rlen,
										Areadnames.begin() + O[bread].nameoff,
										context.fragment,
										secondary,
										supplementary,
										*Pbamheader,
										aux_a,
										aux_e
									);
								}
								else
								{
									context.converter.convert(
										P.first,
										reinterpret_cast<char const *>(refdb.bases) + refdb.reads[aread].boff,
										refdb.reads[aread].rlen,
										context.ARC.begin()+1 /* skip terminator */,
										readlen,
										Areadnames.begin() + O[bread].nameoff,
										context.fragment,
										secondary,
										supplementary,
										*Pbamheader,
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
			}

			for ( uint64_t zr = 1; zr < ointv.size(); ++zr )
			{
				uint64_t const tid = zr-1;
				LASToBamContext & context = *(Acontexts[tid]);
				libmaus2::bambam::parallel::FragmentAlignmentBufferFragment & fragment = context.fragment;

				uint8_t const * bufferstart = fragment.pa;
				uint8_t const * bufferend = fragment.pc;
				char const * cbufferstart = reinterpret_cast<char const *>(bufferstart);
				uint64_t const len = bufferend - bufferstart;
				bgzfout->write(cbufferstart,len);
			}

			if ( OVLdata.size() )
				prevmark = libmaus2::dazzler::align::OverlapData::getBRead(OVLdata.getData(OVLdata.size()-1).first);

			if ( verbose > 1 )
				std::cerr << "[V] processed " << OVLdata.size() << " alignments " << libmaus2::util::MemUsage() << std::endl;
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

	bgzfout->flush();
	bgzfout->addEOFBlock();
	bgzfout.reset();

	Free_Align_Spec(aspec);


	return EXIT_SUCCESS;
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
	ostr << " -i: maximum number of read bases per input block (default 64m)\n";
	ostr << " -M: damapper memory limit\n";
	ostr << " --sasamplingrate: SA sampling rate (default 32)\n";
	ostr << " --bwtconstrmem: memory used to construct BWT (default 3/4 of machine's memory)\n";
	ostr << " -T: prefix for temporary files used during index construction\n";
	ostr << " -Q: file name of index\n";

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
