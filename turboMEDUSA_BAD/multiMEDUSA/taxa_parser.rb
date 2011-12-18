param = 'MedusaPARAM.txt'
taxa = 'MedusaTAXA.txt'
taxact = 0

class TaxaData
  attr_reader :tree, :shift, :taxa, :params

  def initialize(tree, shift, taxa, params)
    @tree = tree
    @shift = shift
    @taxa = taxa
    @params = params
  end

  def to_s
    "#{tree}\t#{shift}\t#{params.join("\t")}"
  end
end

by_taxa = Hash.new

File.open(param) do |pfile|
  File.open(taxa) do |tfile|
    File.open('all_shift_1.txt', 'w') do |s1|
      while true do
        begin
          ps = pfile.readline.split(/\s+/)
          ts = tfile.readline.split(/\s+/, 4)
          td = TaxaData.new(ps.shift, ps.shift, ts[3], ps)

          if td.shift.to_i == 1
            s1.puts(td)
          else
            tx = td.taxa
            by_taxa[tx] ||= Array.new
            by_taxa[tx].push(td)
          end
        rescue EOFError

          by_taxa.each do |k,v|
            File.open("taxa_#{taxact.to_s}", 'w') do |f|
              f.puts(k)
              v.each {|d| f.puts(d)}
            end
            taxact += 1
          end

          exit
        end
      end
    end
  end
end
