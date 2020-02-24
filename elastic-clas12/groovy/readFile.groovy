#!/opt/coatjava/6.3.1/bin/run-groovy

//example code to read a file, check if it exists, and fill a map

p_corr =[:]
for (i in -1..8) {
    sector_corr = []

    def file = new File("p_ele_dp_ele_from_angles_FTOF_1_thetabin${i}_isr.txt")
    if( file.exists() ){
	file.eachLine {  
            line -> println "line : $line"; 
	    def arr = line.tokenize(' ')
	    println("${arr[0]} and ${arr[1]}")
	    sector_corr.add(arr[1])
	}
	p_corr.put(i,sector_corr)
    }    
}

print('done ')
println(p_corr)
