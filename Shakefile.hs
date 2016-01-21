import System.Directory(getCurrentDirectory)
import Control.Monad(liftIO)
import Development.Shake.Command(cmd) 
import Development.Shake
import Development.Shake.FilePath

BIN = "bin" 
ILLUMINA_URL = "http://www.niehs.nih.gov/research/resources/assets/docs/artbinchocolatecherrycake031915linux64tgz.tgz"
PBSIM_URL = "https://pbsim.googlecode.com/files/pbsim-1.0.3-Linux-amd64.tar.gz"

wgetAndUnTar url = do 
  _ <- cmd "wget" url
  let tar = reverse $ takeWhile (/= '/') $ reverse url
  cmd "tar" "-xf" tar

downloadAndLink url binPath linkPath = do 
   _ <- wgetAndUnTar url
   pwd <- liftIO getCurrentDirectory
   cmd "ln" "-s" pwd </> binPath linkPath

BIN </> "art_illumina" %> \out -> do
  downloadAndLink ILLUMINA_URL "art_bin_ChocolateCherryCake/art_illumina" out

BIN </> "pbsim" %> \out -> do
   downloadAndLink PBSIM_URL "pbsim-1.0.3-Linux-amd64/Linux-amd64/bin/pbsim" out

["paired_end_com1.aln" ,"paired_end_com1.fq"
,"paired_end_com2.aln" ,"paired_end_com2.fq"
,"paired_end_com.sam"] &%> \_, _, _, _, _ -> do

  need [BIN </> "art_illumina"]
  cmd "art_illumina" ["-i", "art_bin_ChocolateCherryCake/examples/testSeq.fa", "-o", "./paired_end_com", "-ss", "HS25", "-l", "150", "-f", "10", "-p", "-m", "500", "-s", "10", "-sam"]

["sd_0001.fastq", "sd_0001.maf", "sd_0001.ref"] &%> \_,_,_ -> do

  need [BIN </> "pbsim"]
  cmd "pbsim" ["--data-type", "CLR", "--depth", "20", "--model_qc", "data/model_qc_clr", "sample/sample.fasta"]
