package se.lth.immun

import java.io.File
import scala.io.Source
import scala.collection.mutable.ArrayBuffer

object XLMasterReader {

	
	val COMMON_HEADER = Array(
			"mass", 
			"z", 
			"mz", 
			"xlink"
		)
	
	val TAXLINK_HEADER = Array(
			"lightMonoMz",
			"heavyMonoMz",
			"rtApexLight",
			"rtApexHeavy",
			"massLight",
			"massHeavy",
			"intLight",
			"intHeavy",
			"intScore",
			"rtScore",
			"massScore",
			"scoreSum"
		)
		
	val KOJAK_PERC_HEADER = Array(
			"score",
			"q-value",
			"posteriorErrorProb",
			"psmID",
			"kojakPPM",
			"kojakNormRank",
			"proteinIDs"			
		)
		
		
	case class Row(
			common:Common,
			taxlink:Taxlink,
			kojak:Kojak
		) {
		def isComplete = 
			taxlink.hasResult && kojak.hasResult
	}
	case class Common(
			mass:Double,
			z:Int,
			mz:Double,
			xlink:String
		)
	case class Taxlink(
			lightMonoMz:Double,
			heavyMonoMz:Double,
			rtApexLight:Double,
			rtApexHeavy:Double,
			massLight:Double,
			massHeavy:Double,
			intLight:Double,
			intHeavy:Double,
			intScore:Double,
			rtScore:Double,
			massScore:Double,
			scoreSum:Double
		) {
		def hasResult =
			heavyMonoMz > 0
	}
	case class Kojak(
			score:Double,
			qValue:Double,
			posteriorErrorProb:Double,
			psmID:String,
			kojakPPM:Double,
			kojakNormRank:Int,
			proteinIDs:Seq[String]
		) {
		def hasResult =
			qValue >= 0
	}
	
	def read(path:String):ArrayBuffer[Row] = {
		var headerParsed = false
		val f = new File(path)
		
		var iMASS = -1
		var iZ = -1
		var iMZ = -1
		var iXLINK = -1
		var iLIGHT_MONO_MZ = -1
		var iHEAVY_MONO_MZ = -1
		var iRT_APEX_LIGHT = -1
		var iRT_APEX_HEAVY = -1
		var iMASS_LIGHT = -1
		var iMASS_HEAVY = -1
		var iINT_LIGHT = -1
		var iINT_HEAVY = -1
		var iINT_SCORE = -1
		var iRT_SCORE = -1
		var iMASS_SCORE = -1
		var iSCORE_SUM = -1
		var iSCORE = -1
		var iQ_VALUE = -1
		var iPOST_ERR_PROB = -1
		var iPSM_ID = -1
		var iKOJAK_PPM = -1
		var iKOJAK_NORM_RANK = -1
		var iPROTEIN_IDS = -1
		
		val a = new ArrayBuffer[Row]
		for (line <- Source.fromFile(f).getLines) {
			val cols = line.split("\t")map(_.trim)
			if (!headerParsed) {
				iMASS = cols.indexOf("mass")
				iZ = cols.indexOf("z")
				iMZ = cols.indexOf("mz")
				iXLINK = cols.indexOf("xlink")
				iLIGHT_MONO_MZ = cols.indexOf("lightMonoMz")
				iHEAVY_MONO_MZ = cols.indexOf("heavyMonoMz")
				iRT_APEX_LIGHT = cols.indexOf("rtApexLight")
				iRT_APEX_HEAVY = cols.indexOf("rtApexHeavy")
				iMASS_LIGHT = cols.indexOf("massLight")
				iMASS_HEAVY = cols.indexOf("massHeavy")
				iINT_LIGHT = cols.indexOf("intLight")
				iINT_HEAVY = cols.indexOf("intHeavy")
				iINT_SCORE = cols.indexOf("intScore")
				iRT_SCORE = cols.indexOf("rtScore")
				iMASS_SCORE = cols.indexOf("massScore")
				iSCORE_SUM = cols.indexOf("scoreSum")
				iSCORE = cols.indexOf("score")
				iQ_VALUE = cols.indexOf("q-value")
				iPOST_ERR_PROB = cols.indexOf("posteriorErrorProb")
				iPSM_ID = cols.indexOf("psmID")
				iKOJAK_PPM = cols.indexOf("kojakPPM")
				iKOJAK_NORM_RANK = cols.indexOf("kojakNormRank")
				iPROTEIN_IDS = cols.indexOf("proteinIDs")
				headerParsed = true
			} else {
				a += Row(
						Common(
							cols(iMASS).toDouble,
							cols(iZ).toInt,
							cols(iMZ).toDouble,
							cols(iXLINK)
						),
						Taxlink(
							cols(iLIGHT_MONO_MZ).toDouble,
							cols(iHEAVY_MONO_MZ).toDouble,
							cols(iRT_APEX_LIGHT).toDouble,
							cols(iRT_APEX_HEAVY).toDouble,
							cols(iMASS_LIGHT).toDouble,
							cols(iMASS_HEAVY).toDouble,
							cols(iINT_LIGHT).toDouble,
							cols(iINT_HEAVY).toDouble,
							cols(iINT_SCORE).toDouble,
							cols(iRT_SCORE).toDouble,
							cols(iMASS_SCORE).toDouble,
							cols(iSCORE_SUM).toDouble
						),
						Kojak(
							cols(iSCORE).toDouble,
							cols(iQ_VALUE).toDouble,
							cols(iPOST_ERR_PROB).toDouble,
							cols(iPSM_ID),
							cols(iKOJAK_PPM).toDouble,
							cols(iKOJAK_NORM_RANK).toInt,
							if (iPROTEIN_IDS >= cols.length) Nil else cols(iPROTEIN_IDS).split("|")
						)
					)
			}
		}
		
		a
	}
}