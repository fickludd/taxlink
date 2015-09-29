package se.lth.immun

import java.io.File
import scala.io.Source
import scala.collection.mutable.ArrayBuffer

object TaxlinkOutput {
	
	val NO_PAIR = TaxlinkPair(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, 0, "No feature found")
	case class TaxlinkPair(
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
			scoreSum:Double,
			z:Int,
			xlink:String
		) {
		def avgRtInMinutes = (rtApexHeavy + rtApexLight) / 120
	}

	def read(path:String):ArrayBuffer[TaxlinkPair] = {
		var headerParsed = false
		val f = new File(path)
		
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
		var iZ = -1
		var iXLINK = -1
		
		val a = new ArrayBuffer[TaxlinkPair]
		for (line <- Source.fromFile(f).getLines) {
			val cols = line.split("\t")map(_.trim)
			if (!headerParsed) {
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
				iZ = cols.indexOf("z")
				iXLINK = cols.indexOf("xlink")
				headerParsed = true
			} else {
				a += TaxlinkPair(
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
					cols(iSCORE_SUM).toDouble,
					cols(iZ).toInt,
					cols(iXLINK)
					
				)
			}
		}
		a
	}
}