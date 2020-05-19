$date = Get-Date

$dir =  "D:\Drive\Projetos\Doutorado2C\ValveEstimation\ValveEstimation\x64\Release\"
$out_dir = "D:\Drive\Projetos\Doutorado2C\test_data\"
$out_file = $out_dir + "log_" + $date.Day + "_" + $date.Month + "_" + $date.Year + "-" + $date.Hour + "_" + $date.Minute + "_" + $date.Second + ".txt"

Write-Host ("Logfile location: " + $out_file)

$params = " --directory=" + $out_dir + " --type=estimation --excitation=sinusoidal --models=sgms-gms -k -n 5"

&($dir + "ValveEstimation.exe") $params.Split(" ") | Out-File ($out_file)