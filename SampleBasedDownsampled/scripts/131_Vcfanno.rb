#!/usr/bin/env ruby
require 'optparse'
require 'fileutils'

class Klass
  WORKDIR = File.expand_path("..", __dir__)
  OUTD    = File.join(WORKDIR, "data")
  LOGD    = File.join(WORKDIR, "logs")
  VCFANNO = "/xcatopt/vcfanno-v0.3.5/vcfanno_linux64"
  CONF = File.join(WORKDIR, "Configs", "hg38.toml")
  CORE    = 32
  TEMPSTUFF = File.join(OUTD, "__TEMP__")
  BCFTOOLS = "/xcatopt/bcftools-1.19/bcftools"

  
  attr_reader :opts
  
  def initialize(opts)
    @opts = opts
  end
  
  def run
    FileUtils.mkdir_p(LOGD)
    FileUtils.mkdir_p(TEMPSTUFF)
    Dir["#{OUTD}/pgds1kgp.bisnv.neighbor?.vcf.gz"].sort.each do |vcf|
    puts vcf
      out = File.basename(vcf).sub(/\.(vcf)\.gz$/, ".annot.vcf.gz")
      log = File.join(LOGD, "#{out}.log")
      next if (File.exist?(out) && !opts[:f])
      cmd = Array.new
      cmd << VCFANNO
      cmd << "-p #{CORE}"
      cmd << CONF
      cmd << vcf
      cmd << "2> #{log}"
      cmd << "|"
      cmd << BCFTOOLS
      cmd << "view"
      cmd << "-Oz"
      # write output directly into TEMPSTUFF so final file stays there
      temp_out = File.join(TEMPSTUFF, out)
      cmd << "-o #{temp_out}"
      warn    cmd.join(" ")
      system  cmd.join(" ") unless opts[:n]
    end
  end
end

if $0 == __FILE__
  opts = ARGV.getopts("hnf")
  Klass.new(opts.map{|k,v|[k.to_sym, v]}.to_h).run
end
