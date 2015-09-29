package se.lth.immun

import java.io.File
import java.io.FileReader
import java.io.FileInputStream
import java.io.BufferedReader
import java.io.InputStreamReader
import java.util.zip.GZIPInputStream

import se.lth.immun.xml.XmlReader
import se.lth.immun.mzml._
import se.lth.immun.mzml.ghost._

object MzMLBasedRTMap {

	def getReader(path:String):XmlReader = {
		val f = new File(path)
		if (path.toLowerCase.endsWith(".mzml.gz"))
			new XmlReader(new BufferedReader(new InputStreamReader(
							new GZIPInputStream(new FileInputStream(f)))))
		else if (path.toLowerCase.endsWith(".mzml"))
			new XmlReader(new BufferedReader(new FileReader(f)))
		else
			throw new Exception("Unknown file format '%s'".format(path))
	}
	
	def read(path:String, params:TaxPerKojakParams):Int => Double = {
		val r = getReader(path)
		var rts = Array[Double]()
		var iSpec = 0
		val t0 = System.currentTimeMillis
		println(" reading retention times from mzML...")
		val dh = new MzMLDataHandlers(
				nSpec => {
					rts = new Array[Double](nSpec)
					println("  n spectra: " + nSpec)
					if (params.verbose)
						println(" scanNum  msLevel  rt    time elapsed")
				},
				spec => {
					val gs = GhostSpectrum.fromSpectrum(spec)
					rts(iSpec) = gs.scanStartTime
					if (iSpec % 500 == 0 && params.verbose) {
						println("%8d %8d %10.3f %ds".format(iSpec, gs.msLevel, gs.scanStartTime, (System.currentTimeMillis - t0) / 1000))
					}
					iSpec += 1
				},
				nChrom => {},
				chrom => {})
		
		MzML.fromFile(r, dh)
		println(" rts read.")
		rts
	}
}