package se.lth.immun

import java.util.Properties
import java.io.File
import se.jt.CLIApp

import scala.collection.mutable.ArrayBuffer
import scala.collection.mutable.HashMap
import scala.collection.mutable.HashSet

object TaxPerKojak extends CLIApp {
	
	import TaxlinkOutput.TaxlinkPair
	import KojakOutput.KojakRow
	import KojakPercOutput.PercRow
	import MasterWriter.Row
	
	case class KojakPercRow(kojakRow:Option[KojakRow], percRow:Option[PercRow])
	
	def main(args:Array[String]):Unit = {
		var properties = new Properties
	    properties.load(this.getClass.getResourceAsStream("/pom.properties"))
		val appName 	= properties.getProperty("pom.artifactId")
		val appVersion 	= properties.getProperty("pom.version")
    
		implicit val params = new TaxPerKojakParams
		
		failOnError(parseArgs(appName, appVersion, args, params, List(), None))
		
		if (params.percOpt.isEmpty && params.taxOpt.isEmpty && params.kojakOpt.isEmpty) {
			print(usage(appName, appVersion, args, params, List(), None))
			println("Need to give at least one of 'taxlink' or 'kojak' and 'perc'!")
			System.exit(1)
		}
		
		println("   "+appName+" "+appVersion)
		for (path <- params.taxOpt)
			println("      taxlink file: " + path)
		for (path <- params.kojakOpt)
			println("        kojak file: " + path)
		for (path <- params.percOpt)
			println("   percolator file: " + path)
		val taxlinkPairs 	= params.taxOpt.map(TaxlinkOutput.read)
		val percRows 		= params.percOpt.map(KojakPercOutput.read)
		val kojakRows 		= params.kojakOpt.map(KojakOutput.read)

		
		
		val outRows:Seq[Row] = 
			(taxlinkPairs, kojakRows, percRows) match {
				case (None, None, None) =>
					throw new Exception("This should have been caught earlier!")
					Nil
					
				case (Some(tps), None, None) =>
					tps.map(tp => Row(Some(tp), None, None))
		
				case (_, Some(krows), None) =>
					failOnError(Array("This is stupid... run percolator for statistics!"))
					Nil

				case (_, None, Some(prows)) =>
					failOnError(Array("This is not possible because of charge state ambiguity... include kojak output!"))
					Nil
					
				case (None, Some(krows), Some(prows)) =>
					joinKojakFiles(krows, prows).map(kp => Row(None, Some(kp), None))
				
				case (Some(tps), Some(krows), Some(prows)) =>
					
					val rtMapOpt = 
						if (krows.forall(_.rt.nonEmpty)) None
						else 
							params.mzMLOpt match {
								case Some(path) => 
									println("         mzML file: " + path)
									Some(MzMLBasedRTMap.read(path, params))
								case None => 
									println("WARNING: without the mzML kojak and taxlink result can't be matched based on RT. This could be dangerous.")
									None
							}
					
					val kojakPercRows = joinKojakFiles(krows, prows)
					val pmap = new HashMap[String, ArrayBuffer[(KojakRow, PercRow)]]
					for (kp <- kojakPercRows) {
						val xl = kp._2.peptide
						if (!pmap.contains(xl))
							pmap += xl -> new ArrayBuffer
						pmap(xl) += kp
					}
					val tpBased = for (tp <- tps) yield {
						val kpOpt = pmap.get(tp.xlink).flatMap(kps =>
								kps.find(kp => taxlinkKojakMatch(tp, kp._1, params, rtMapOpt))
							)
						val rtOpt =
							for {
								kp <- kpOpt
								rtMap <- rtMapOpt
							} yield rtMap(kp._1.scanNum)
						Row(Some(tp), kpOpt, rtOpt)
					}
					val matchedKPs = tpBased.flatMap(_.kojakPercRow).toSet
					val kpBased = 
						for {
							kps <- pmap.values
							kp <- kps
							if !matchedKPs.contains(kp)
						} yield Row(None, Some(kp), rtMapOpt.map(rtMap => rtMap(kp._1.scanNum)))
					
					tpBased ++ kpBased
		}
		
		val (outDir, outName) = params.outBase
		val qvalFiltered = outRows.filter(_.kojakPercRow.map(_._2.qValue).getOrElse(-1.0) < params.qValue)
		
		println("       %d out rows of which ".format(qvalFiltered.length))
		println("            complete: " + qvalFiltered.count(r => r.taxlinkPair.nonEmpty && r.kojakPercRow.nonEmpty))
		println("        taxlink only: " + qvalFiltered.count(_.taxlinkPair.nonEmpty))
		println("          kojak only: " + qvalFiltered.count(_.kojakPercRow.nonEmpty))
		
		MasterWriter.write(new File(outDir, outName+".xl.tsv"), qvalFiltered, params)
	}
	
	
	def taxlinkKojakMatch(tp:TaxlinkPair, k:KojakRow, params:TaxPerKojakParams, rtMapOpt:Option[Int => Double]):Boolean = 
		tp.z == k.charge && 
		withinPPM(tp.massLight, k.mass, params.matchPPM) &&
		rtMapOpt.map(rtMap => math.abs(rtMap(k.scanNum) - tp.avgRtInMinutes) < params.matchRT).getOrElse(true)
	
	
	def withinPPM(m1:Double, m2:Double, ppm:Double) = 
		math.abs(m1 - m2) * 1e6 / m2 < ppm

	
	def joinKojakFiles(
			kojakRows:Seq[KojakRow], 
			percRows:Seq[PercRow]
	):Seq[(KojakRow, PercRow)] = {
		val kmap = new HashMap[String, KojakRow]
		kmap ++= kojakRows.map(krow => krow.specId -> krow)
		percRows.map(prow => (kmap(prow.PSMId), prow))
	}
}