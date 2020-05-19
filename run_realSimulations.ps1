

$dir =  "D:\Drive\Projetos\Doutorado2C\ValveEstimation\ValveEstimation\x64\Release\"
$out_dir = "D:\Drive\Projetos\Doutorado2C\test_data\"

$valves = ("graphite")
$excitations = ("sinusoidal", "aleatory", "ramps")

$params = " --directory=" + $out_dir + " --type=real-simulation --models=he-choudhury-kano-karnopp-lugre-gms-sgms-gms1"
foreach ($valve in $valves){
	foreach ($excitation in $excitations){
		$params += " --excitation=" + $excitation + " --valve=" + $valve + " --load=D:\Drive\Projetos\Doutorado2C\test_data\" + $valve + "_data_ol_" + $excitation + ".csv"
		&($dir + "ValveEstimation.exe") $params.Split(" ")
	}
}



# $dir =  "D:\Drive\Projetos\Doutorado2C\ValveEstimation\ValveEstimation\x64\Release\"
# $out_dir = "D:\Drive\Projetos\Doutorado2C\test_data\"

# $valves = ("graphite", "teflon")
# $excitations = ("sinusoidal")

# $params = " --directory=" + $out_dir + " --type=real-simulation --models=he-choudhury-kano-karnopp-lugre-gms-sgms-gms1"
# foreach ($valve in $valves){
	# foreach ($excitation in $excitations){
		# $params += " --excitation=" + $excitation + " --valve=" + $valve + " --load=D:\Drive\Projetos\Doutorado2C\test_data\" + $valve + "_data_ol_" + "ramps" + ".csv"
		# &($dir + "ValveEstimation.exe") $params.Split(" ")
	# }
# }