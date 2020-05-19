

$dir =  "D:\Drive\Projetos\Doutorado2C\ValveEstimation\ValveEstimation\x64\Release\"
$out_dir = "D:\Drive\Projetos\Doutorado2C\test_data\"

$valves = ("graphite")
$excitations = ("sinusoidal")

$params = " --directory=" + $out_dir + " --type=real-estimation --models=he -k"
foreach ($valve in $valves){
	foreach ($excitation in $excitations){
		$date = Get-Date
		$out_file = $out_dir + "log_" + $date.Day + "_" + $date.Month + "_" + $date.Year + "-" + $date.Hour + "_" + $date.Minute + "_" + $date.Second + ".txt"
		$params += " --excitation=" + $excitation + " --valve=" + $valve + " --load=D:\Drive\Projetos\Doutorado2C\test_data\" + $valve + "_data_ol_" + $excitation + ".csv" + " -n 5"
		&($dir + "ValveEstimation.exe") $params.Split(" ") | Out-File $out_file
	}
}