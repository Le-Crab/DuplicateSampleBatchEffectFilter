#!/usr/bin/env ruby
require 'optparse'
require 'bio'
require 'bio/twobit'
require 'statistics2'

class Klass
  WORKDIR = File.expand_path("..", __dir__)
  BCFTOOLS = "/xcatopt/bcftools-1.18/bcftools"
  NEIGHBOR_SIZE = 2
  VCFIN  = File.join(WORKDIR, "data", "pgds1kgp.bisnv.bcf")
  VCFOUT = VCFIN.sub(/\.bcf\z/, ".neighbor#{NEIGHBOR_SIZE}.vcf.gz")
  #  TWOBIT = "/peta/fserv/share/Genomes/human_hg38_GRCh38/hg38.2bit"
  TWOBIT = "/peta/fserv/share/Genomes/human_hg38_GRCh38/hg38.canonical-EBV-JRGv2.2bit"

  NUMADJ = 2

  INFO_IX   = 7
  SAMPLE_IX = 9
  PHRED_MAX = Float::MIN_10_EXP * -10

  OMIT = "__OMIT_ME__"


  NUCLEOTIDES = %w(A C G T)

  ADDITIONAL_HEADERS =
    ['##INFO=<ID=GINeighbors,Number=1,Type=String,Description="Substituion and adjacent nucleotides">',
     '##INFO=<ID=GIConcordant,Number=0,Type=Flag,Description="Concordant among paired samples">',
     '##INFO=<ID=GIDiscordant,Number=0,Type=Flag,Description="Discordant among paired samples">',
    ]

  PAIRED = [
    # https://www.internationalgenome.org/data-portal/data-collection/platinum-pedigree
    # Platinum    1000genomes 
    ["ERR194158", "NA12889"],
    ["ERR194159", "NA12890"],
    ["ERR194160", "NA12891"],
    ["ERR194161", "NA12892"],
    ["ERR194146", "NA12877"], # child of trio #1
    ["ERR194147", "NA12878"], # child of trio #2
  ]
  
  attr_reader :opts, :paired_ixs
  
  def initialize(opts)
    @opts = opts
  end


  def twobit
    @twobit ||= Bio::TwoBit.open(TWOBIT)
    @twobit
  end

  def complement_neighbor(seq0)
    seq = seq0.gsub(/:/,'')
    ss = seq.size
    fail "Not expecpted length of the neighbors-code" unless ss.size != (NEIGHBOR_SIZE * 2) + 2 
    subst_idx = (ss / 2 ) - 1
    unless @compnei
      @compnei = Hash.new
      perms =
        NUCLEOTIDES.repeated_permutation(ss).map(&:join).reject{|x|x[subst_idx] == x[subst_idx + 1]}
      perms.reverse_each do |perm|
        rcomp = Bio::Sequence::NA.new(perm).reverse_complement.upcase
        (perm, rcomp = rcomp, perm) if perm > rcomp
        @compnei[perm]  = perm
        @compnei[rcomp] = perm
      end
      # When (NEIGHBOR_SIZE == 2)
      # 3072 patterns -> 1536
      # AAACAA > AAACAA
      # AAACAC > AAACAC
      # TTTGCG > CGCAAA
      # TTTGCT > AGCAAA
    end
    res = @compnei[seq]
    res = "N" * ((NEIGHBOR_SIZE * 2) + 2) if res.nil?
    left  = res[0..(NEIGHBOR_SIZE - 1)]
    subst = res[(NEIGHBOR_SIZE) .. (NEIGHBOR_SIZE + 1 )]
    right = res[(NEIGHBOR_SIZE + 2) .. -1]
    "#{left}:#{subst}:#{right}"
  end

  def neighbors(cols, comp:)
    chrom = cols[0]
    pos_st = Integer(cols[1]) - 1 # 0-based [st, ed)
    pos_ed = Integer(cols[1])     # 0-based [st, ed)
    adj_bef = [chrom, (pos_st - NUMADJ), pos_st]
    #adj_aft = [chrom, (pos_ed - 1), (pos_ed - 1 + NUMADJ)]
    adj_aft = [chrom, pos_ed, (pos_ed + NUMADJ)]
    adj_bef_nuc = twobit.sequence(*adj_bef)
    adj_aft_nuc = twobit.sequence(*adj_aft)
    # "GINeighbors=#{adj_bef_nuc}:#{cols[3]}:#{cols[4]}:#{adj_aft_nuc}"
    encode = "#{adj_bef_nuc}:#{cols[3]}:#{cols[4]}:#{adj_aft_nuc}"
    cneighbor = complement_neighbor(encode)
    if cneighbor.include?("N")
       OMIT
    else
      if comp
        # return "GICompNeighbors=#{complement_neighbor(encode)}"
      else
        return "GINeighbors=#{encode}"
      end
    end
  end
  

  def find_paired_indexs(titlecols)
    @paired_ixs ||= Array.new
    samples = titlecols[SAMPLE_IX .. -1]
    PAIRED.each do |pair|
      @paired_ixs <<
        [ samples.index{|x|x==pair[0]}, samples.index{|x|x==pair[1]} ]
    end
  end

  def concdisc(cols)
    gts = cols[SAMPLE_IX .. -1].map{|s|s.split(":")[0]}
    return OMIT if paired_ixs.any?{|(ix1, ix2)| gts[ix1].include?(".") || gts[ix2].include?(".")} 
    if paired_ixs.all?{|(ix1, ix2)| gts[ix1] == gts[ix2]}
      "GIConcordant"
    else
      "GIDiscordant"
    end
  end
  
  def run
    IO.popen("#{BCFTOOLS} view -Oz -o #{VCFOUT} -", "w") do |fout|
      IO.popen("#{BCFTOOLS} view #{VCFIN}", "r") do |fin|
        fin.each_line do |row|
          if row.start_with?("##")
            fout.puts row
            next
          end
          if row.start_with?("#CHROM")
            fout.puts ADDITIONAL_HEADERS.join("\n")
            fout.puts row
            find_paired_indexs(row.chomp.split("\t"))  
            next
          end
          cols = row.chomp.split("\t")
          next if (cols[3].length > 1) || (cols[4].length > 1) # skip MNV like ATC>TGC
          cols[INFO_IX] << ";#{neighbors(cols, comp:false)}"
          cols[INFO_IX] << ";#{neighbors(cols, comp:true)}"
          cols[INFO_IX] << ";#{concdisc(cols)}"
          next if cols[INFO_IX].include?(OMIT)
          fout.puts cols.join("\t")
        end
      end
    end
  end
end

if $0 == __FILE__
  opts = ARGV.getopts("h")
  Klass.new(opts.map{|k,v|[k.to_sym, v]}.to_h).run
end
