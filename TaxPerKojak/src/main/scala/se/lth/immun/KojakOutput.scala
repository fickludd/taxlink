package se.lth.immun

import java.io.File
import scala.io.Source
import scala.collection.mutable.ArrayBuffer

object KojakOutput {

	val NO_ROW = KojakRow("", "", -1, None, -1, -1, -1, -1, 0, -1, -1, 0,0,0,"",Nil)
	
	case class KojakRow(
			specId:String,
			label:String,
			scanNum:Int,
			rt:Option[Double],
			score:Double,
			dScore:Double,
			normRank:Int,
			ppScoreDiff:Double,
			charge:Int,
			mass:Double,
			ppm:Double,
			lenShort:Int,
			lenLong:Int,
			lenSum:Int,
			peptide:String,
			proteins:Seq[String]
			)
	
	val SINGLE_LOOP_HEADER = "SpecId	Label	scannr	Score	dScore	Charge	Mass	PPM	Len	Peptide	Proteins"
	val INTRA_INTER_HEADER = "SpecId	Label	scannr	Score	dScore	NormRank	PPScoreDiff	Charge	Mass	PPM	LenShort	LenLong	LenSum	Peptide	Proteins"
	
	def read(path:String):ArrayBuffer[KojakRow] = {
		var headerParsed = false
		val f = new File(path)
		
		var iSPEC_ID = -1
		var iLABEL = -1
		var iSCAN_NUM = -1
		var iSCORE = -1
		var iD_SCORE = -1
		var iNORM_RANK = -1
		var iPP_SCORE_DIFF = -1
		var iCHARGE = -1
		var iMASS = -1
		var iPPM = -1
		var iLEN_SHORT = -1
		var iLEN_LONG = -1
		var iLEN_SUM = -1
		var iPEPTIDE = -1
		var iPROTEINS = -1
		
		val a = new ArrayBuffer[KojakRow]
		
		for (line <- Source.fromFile(f).getLines) {
			val cols = line.split("\t")map(_.trim)
			if (!headerParsed) {
				if (line == SINGLE_LOOP_HEADER) {
					iSPEC_ID = 0
					iLABEL = 1
					iSCAN_NUM = 2
					iSCORE = 3
					iD_SCORE = 4
					iCHARGE = 5
					iMASS = 6
					iPPM = 7
					iLEN_SHORT = 8
					iLEN_LONG = 9
					iLEN_SUM = 10
					iPEPTIDE = 11
					iPROTEINS = 12
				} else { // presumably INTRA_INTER, but the general parser should work
					iSPEC_ID = cols.indexOf("SpecId")
					iLABEL = cols.indexOf("Label")
					iSCAN_NUM = cols.indexOf("ScanNum")
					iSCORE = cols.indexOf("Score")
					iD_SCORE = cols.indexOf("dScore")
					iNORM_RANK = cols.indexOf("NormRank")
					iPP_SCORE_DIFF = cols.indexOf("PPScoreDiff")
					iCHARGE = cols.indexOf("Charge")
					iMASS = cols.indexOf("Mass")
					iPPM = cols.indexOf("PPM")
					iLEN_SHORT = cols.indexOf("LenShort")
					iLEN_LONG = cols.indexOf("LenLong")
					iLEN_SUM = cols.indexOf("LenSum")
					iPEPTIDE = cols.indexOf("Peptide")
					iPROTEINS = cols.indexOf("Proteins")
				}
				headerParsed = true
			} else {
				val specId = cols(iSPEC_ID)
				val scanNum = cols(iSCAN_NUM).toInt
				a += KojakRow(
					specId,
					cols(iLABEL),
					scanNum, 
					parseSpecIdRt(specId),
					cols(iSCORE).toDouble,
					cols(iD_SCORE).toDouble,
					if (iNORM_RANK < 0) -1 else cols(iNORM_RANK).toInt,
					if (iPP_SCORE_DIFF < 0) -1 else cols(iPP_SCORE_DIFF).toDouble,
					cols(iCHARGE).toInt,
					cols(iMASS).toDouble,
					cols(iPPM).toDouble,
					cols(iLEN_SHORT).toInt,
					cols(iLEN_LONG).toInt,
					cols(iLEN_SUM).toInt,
					cols(iPEPTIDE),
					cols.drop(iPROTEINS)
				)
			}
		}
		a
	}
	
	def parseSpecIdRt(specId:String) = {
		if (specId.contains("."))
			Some(specId.split("-")(2).toDouble)
		else None
	}
}