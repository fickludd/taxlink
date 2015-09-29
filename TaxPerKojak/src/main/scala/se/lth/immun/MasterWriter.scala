package se.lth.immun

import java.io.File
import java.io.FileWriter
import java.io.BufferedWriter
import se.lth.immun.chem.Constants

object MasterWriter {

	import TaxlinkOutput.TaxlinkPair
	import KojakOutput.KojakRow
	import KojakPercOutput.PercRow
	
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
			"scanNum",
			"rt",
			"score",
			"q-value",
			"posteriorErrorProb",
			"psmID",
			"kojakPPM",
			"kojakNormRank",
			"proteinIDs"			
		)
		
	case class Row(
			taxlinkPair:Option[TaxlinkPair], 
			kojakPercRow:Option[(KojakRow, PercRow)],
			mzMLRt:Option[Double]
		)
	
	def write(f:File, rows:Seq[Row], params:TaxPerKojakParams) = {
		val w = new BufferedWriter(new FileWriter(f))
		
		def writeRow(qoute:Boolean)(a:Any*) = w.write(rowPart(qoute)(a:_*) + "\n")
		def writeRowPart(qoute:Boolean)(a:Any*) = w.write(rowPart(qoute)(a:_*) + params.outSep)
		def rowPart(qoute:Boolean)(a:Any*) = 
			a.map(_ match {
				case s:String => 
					if (qoute) params.outQuote + s + params.outQuote
					else s
				case x => x.toString
			}).mkString(params.outSep)
			
		
			
		writeRowPart(false)(COMMON_HEADER:_*)
		writeRowPart(false)(TAXLINK_HEADER:_*)
		writeRow(false)(KOJAK_PERC_HEADER:_*)
			
		for (Row(tpOpt, kpOpt, mzMLRtOpt) <- rows) {
			val (mass,z,xlink) = 
				tpOpt match {
					case Some(tp) =>
						(tp.massLight, tp.z, tp.xlink)
					case None =>
						val (k, p) = kpOpt.get
						(k.mass, k.charge, k.peptide)
				}
			val mz = (mass + Constants.PROTON_WEIGHT*z) / z
			writeRowPart(false)(mass, z, mz, xlink)
			
			val tp = tpOpt.getOrElse(TaxlinkOutput.NO_PAIR)
			writeRowPart(false)(
					tp.lightMonoMz,
					tp.heavyMonoMz,
					tp.rtApexLight,
					tp.rtApexHeavy,
					tp.massLight,
					tp.massHeavy,
					tp.intLight,
					tp.intHeavy,
					tp.intScore,
					tp.rtScore,
					tp.massScore,
					tp.scoreSum)
					
			val (k, p) = kpOpt.getOrElse((KojakOutput.NO_ROW, KojakPercOutput.NO_ROW))
			writeRow(false)(
					k.scanNum,
					k.rt.orElse(mzMLRtOpt).getOrElse(-1.0),
					p.score,
					p.qValue,
					p.posterior_error_prob,
					p.PSMId,
					k.ppm,
					k.normRank,
					k.proteins.mkString("|")
					)
		}
		
		w.close
	}
}