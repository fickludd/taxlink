package se.lth.immun

class Feature(
		val rt:Double,
		val mz:Double,
		val intensity:Double,
		val charge:Int,
		val width:Double,
		val quality:Double,
		val rt_quality:Double,
		val mz_quality:Double,
		val rt_start:Double,
		val rt_end:Double
) {
	
	override def toString = "%.5f +%d %.1fs".format(mz, charge, rt)
}
