import System.Directory(getCurrentDirectory)
import Control.Monad.IO.Class(liftIO)
import Development.Shake.Command(cmd) 
import Development.Shake
import Development.Shake.FilePath

bin = "bin" 
illuminaUrl = Url "http://www.niehs.nih.gov/research/resources/assets/docs/artbinchocolatecherrycake031915linux64tgz.tgz"
pbsimUrl = Url "https://pbsim.googlecode.com/files/pbsim-1.0.3-Linux-amd64.tar.gz"

illuminaTestOut = ["paired_end_com1.aln" , "paired_end_com1.fq" , "paired_end_com2.aln" , "paired_end_com2.fq" , "paired_end_com.sam"]

newtype Url = Url String

wgetAndUnTar :: Url -> Action () 
wgetAndUnTar (Url url) = do 
  () <- cmd "wget" url
  let tar = reverse $ takeWhile (/= '/') $ reverse url
  cmd "tar" "-xf" tar
  
downloadAndLink :: Url -> FilePath -> FilePath -> Action ()
downloadAndLink url binPath linkPath = do 
   _ <- wgetAndUnTar url
   pwd <- liftIO getCurrentDirectory
   cmd "ln" "-s" (pwd </> binPath) linkPath

main :: IO ()
main = shakeArgs shakeOptions{shakeFiles="_build"} $ do
  want illuminaTestOut
  bin </> "art_illumina" %> \out -> do
    downloadAndLink illuminaUrl "art_bin_ChocolateCherryCake/art_illumina" out
  
  bin </> "pbsim" %> \out -> do
     downloadAndLink pbsimUrl "pbsim-1.0.3-Linux-amd64/Linux-amd64/bin/pbsim" out
  
  illuminaTestOut &%> \_  -> do 
    need [bin </> "art_illumina"]
    cmd (bin </> "art_illumina") ["-i", "art_bin_ChocolateCherryCake/examples/testSeq.fa", "-o", "./paired_end_com", "-ss", "HS25", "-l", "150", "-f", "10", "-p", "-m", "500", "-s", "10", "-sam"]
  
  ["sd_0001.fastq", "sd_0001.maf", "sd_0001.ref"] &%> \_ -> do
  
    need [bin </> "pbsim"]
    cmd "pbsim" ["--data-type", "CLR", "--depth", "20", "--model_qc", "data/model_qc_clr", "sample/sample.fasta"]
