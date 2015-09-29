package se.lth.immun

import java.io.File
import scala.io.Source
import scala.collection.mutable.ArrayBuffer

object KojakPercOutput {

	val NO_ROW = PercRow("", -1, -1, -1, "", Nil)
	
	case class PercRow(
			PSMId:String,
			score:Double,
			qValue:Double,
			posterior_error_prob:Double,
			peptide:String,
			proteinIds:Seq[String]
			)
			
	def read(path:String):ArrayBuffer[PercRow] = {
		var headerParsed = false
		val f = new File(path)
		
		var iPSM_ID = -1
		var iSCORE = -1
		var iQ_VALUE = -1
		var iPOST_ERR_PROB = -1
		var iPEPTIDE = -1
		var iPROTEIN_IDS = -1
		
		val a = new ArrayBuffer[PercRow]
		for (line <- Source.fromFile(f).getLines) {
			val cols = line.split("\t")map(_.trim)
			if (!headerParsed) {
				iPSM_ID = cols.indexOf("PSMId")
				iSCORE = cols.indexOf("score")
				iQ_VALUE = cols.indexOf("q-value")
				iPOST_ERR_PROB = cols.indexOf("posterior_error_prob")
				iPEPTIDE = cols.indexOf("peptide")
				iPROTEIN_IDS = cols.indexOf("proteinIds")
				headerParsed = true
			} else {
				a += PercRow(
					cols(iPSM_ID),
					cols(iSCORE).toDouble,
					cols(iQ_VALUE).toDouble,
					cols(iPOST_ERR_PROB).toDouble,
					cols(iPEPTIDE),
					cols.drop(iPROTEIN_IDS)
				)
			}
			
		}
		a
	}
}