#!/usr/bin/env ruby
require 'optparse'
require 'fileutils'


class Klass
  WORKDIR = File.expand_path("..", __dir__)
  OUTD = File.join(WORKDIR, "data")
  LOGD = File.join(WORKDIR, "logs")
  JAR = "/xcatopt/snpEff-5.1/SnpSift.jar"
  BCFTOOLS = "/xcatopt/bcftools-1.19/bcftools"
  TEMPSTUFF = File.join(OUTD, "__TEMP__")
  attr_reader :opts
  
  def initialize(opts)
    @opts = opts
  end
  
  def run
    FileUtils.mkdir_p(LOGD)
    puts(WORKDIR)
    Dir["#{OUTD}/*.neighbor2.annot.vcf.gz"].sort.each do |vcf|
      out_base = File.basename(vcf).sub(/\.(vcf|bcf)\.gz$/, ".txt")
      out = File.join(OUTD, out_base)
      temp_out = File.join(OUTD, "__TEMP__#{out_base}")
      log = File.join(LOGD, "#{out_base}.log")
      pp out
      next if (File.exist?(out) && !opts[:f])
      samples = `#{BCFTOOLS} query --list-samples #{vcf}`.chomp.split("\n")
      #open("#{out}#{TEMPSUFF}", "w"){|fout|fout.puts("#" + samples.join("\t"))} unless opts[:n]
      cmd = Array.new
      cmd << "java -Xmx128g -jar #{JAR}"
      cmd << "extractFields"
      cmd << '-s ","'
      cmd << '-e "."'
      cmd << vcf
      cmd << '"CHROM"'
      cmd << '"POS"'
      cmd << '"REF"'
      cmd << '"ALT"'
      #cmd << '"GIMeanDepth"'
      #cmd << '"GIMeanDepthAlt"'
      cmd << '"GINeighbors"'
      #cmd << '"GICompNeighbors"'
      #cmd << '"GIStrandBias"'
      #cmd << '"GIConcordant"'
      #cmd << '"GIDiscordant"'
      cmd << '"RepeatMasker"'
      cmd << '"WindowMaskerSdust"'
      cmd << '"SegDup"'
      cmd << '"CpGisland"'
      cmd << '"SimpleRepeat"'
      cmd << '"Mappability"'
      samples.each do |sname|
        cmd << "\"GEN[#{sname}].GT\""
        cmd << "\"GEN[#{sname}].DP\""
        cmd << "\"GEN[#{sname}].AD\""
        cmd << "\"GEN[#{sname}].SB[*]\""
        cmd << "\"GEN[#{sname}].GQ\""
        #cmd << "\"GEN['#{sname}'].PL\""
        #cmd << "\"GEN['#{sname}'].RNC\""
      end
      ###cmd << ">> #{out}#{TEMPSUFF}"
  cmd << "> #{temp_out}"
      cmd << "2> #{log}"
      cmd << "&&"
  cmd << "mv #{temp_out} #{out}"
      warn   cmd.join(" ")
      system cmd.join(" ") unless opts[:n]
    end
  end
end

if $0 == __FILE__
  opts = ARGV.getopts("hnf")
  Klass.new(opts.map{|k,v|[k.to_sym, v]}.to_h).run
end
