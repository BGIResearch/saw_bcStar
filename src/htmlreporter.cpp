#include "htmlreporter.h"
#include <chrono>
#include <memory.h>

extern string command;

HtmlReporter::HtmlReporter(Options* opt){
    mOptions = opt;
}

HtmlReporter::~HtmlReporter(){
}

void HtmlReporter::outputRow(ofstream& ofs, string key, long v) {
    ofs << "<tr><td class='col1'>" + key + "</td><td class='col2'>" + to_string(v) + "</td></tr>\n";
}

void HtmlReporter::outputRow(ofstream& ofs, string key, string v) {
    ofs << "<tr><td class='col1'>" + key + "</td><td class='col2'>" + v + "</td></tr>\n";
}

string HtmlReporter::formatNumber(long number) {
    double num = (double)number;
    string unit[6] = {"", "K", "M", "G", "T", "P"};
    int order = 0;
    while (num > 1000.0) {
        order += 1;
        num /= 1000.0;
    }

    if (order == 0)
        return to_string(number);
    else
        return to_string(num) + " " + unit[order];
}

string HtmlReporter::toThousands(long number){
    std::string numberStr = std::to_string(number);
    int len = numberStr.length();
    for (int index = len-3 ; index > 0; index-=3){
        numberStr.insert(index, ",");
    }
    return numberStr;
}

string HtmlReporter::getPercents(long numerator, long denominator) {
    if(denominator == 0)
        return "0.0";
    else{
        std::ostringstream percentOut;
        percentOut << std::fixed << std::setprecision(2) << (double)numerator * 100.0 / (double)denominator;
        return percentOut.str();
    }
}

void HtmlReporter::reportBarcode(ofstream& ofs, long* barcodeStat, int maxReadsNumber) {
	long* x = new long[maxReadsNumber];
	for (int i = 0; i < maxReadsNumber; i++) {
		x[i] = i +1;
	}
    ofs << resetiosflags(ios::fixed) << setprecision(2);
	ofs << "<div id='barcode_figure'>\n";
	ofs << "<div class='figure' id = 'plot_barcode_reads_number', style='height:400px;'></div>\n";
	ofs << "</div>\n";
	
	ofs << "<div class='sub_section_tips'>This is the statistic of reads number per barcode. <br /> The x axis represents the reads number per barcode while the y axis represent the total reads number corresponding to the barcode with reads number marked by x axis.";
	ofs << "</div>\n";
	ofs << "\n<script type=\"text/javascript\">" << endl;
	string json_str = "var data=[";

	json_str += "{";
	json_str += "x:[" + list2string(x, maxReadsNumber) + "],";
	json_str += "y:[" + list2string(barcodeStat, maxReadsNumber) + "],";
	json_str += "name: 'Reads_Number ',";
	json_str += "type:'bar',";
	json_str += "line:{color:'rgba(225,0,128,1.0)', width:1.5}\n";
	json_str += "}";

	json_str += "];\n";

	json_str += "var layout={title:'barcode reads number ', xaxis:{title:'Reads number per barcode (paired)'}, yaxis:{title:'Reads_number * barcode_number'}};\n";
	json_str += "Plotly.newPlot('plot_barcode_reads_number', data, layout);\n";

	ofs << json_str;
	ofs << "</script>" << endl;

	delete[] x;
}

void HtmlReporter::printReport(string htmlFile, BarcodePositionMap& bpmap){
    ofstream ofs;
    ofs.open(htmlFile, ifstream::out);
    ofs << fixed << setprecision(2);
    
    printHeader(ofs);
    printSummary(ofs, bpmap);
    printFooter(ofs);
    ofs.close();
}

void HtmlReporter::printSummary(ofstream& ofs, BarcodePositionMap& bpmap){
    ofs << endl;
    ofs << "<div style='color:black; background-color: #2696ec;'>\n";
    ofs << "<h1 style='text-align:center; font-family:sans-serif;'>Stereomics barcode QC</h1>\n";
    //ofs << "<h2 style='text-align: left; color: #ffffff; font-family: monospace;'>BGI BIGDATA</h2>\n";
    ofs << "</div>";

    ofs << "<div class='section_title' onclick=showOrHide('SequencingQuality')>Sequencing Quality</div>\n";
    ofs << "<div id='SequencingQuality'>\n";
    ofs << "<div class='subsection_title' onclick=showOrHide('basic')>Basic Information</div>\n";
    ofs << "<div id='basic'>\n";
    ofs << "<table class='summary_table'>\n";
    
    outputRow(ofs, "get unique barcode types:", formatNumber(bpmap.getBarcodeTypes()));
    outputRow(ofs, "barcode reverse complement:", mOptions->rcString);
    outputRow(ofs, "barcode start:", bpmap.barcodeStart);
    outputRow(ofs, "barcode length:", bpmap.barcodeLen);
    outputRow(ofs, "barcode segment:", bpmap.segment);

    ofs << "</table>\n";
    ofs << "</div>\n";
    
    string polymsg = "barcodes with poly";
    for(int i = 0; i< bpmap.inFastqNumber; i++){
        std::string filename = basename(bpmap.inFile[i]);
        ofs << "<div class='subsection_title' onclick=showOrHide('fastqStat')>" + filename + "</div>\n";
        ofs << "<div id='fastqStat'>\n";
        ofs << "<table class='summary_table'>\n";
        outputRow(ofs, "total reads:", formatNumber(bpmap.totalReads[i]));
        outputRow(ofs, "total bases:", formatNumber(bpmap.totalBase[i]));
        outputRow(ofs, "reads without position:", formatNumber(bpmap.readsWithoutPos[i]) + " (" + getPercents(bpmap.readsWithoutPos[i], bpmap.totalReads[i]*bpmap.segment) + "%)");
        outputRow(ofs, "barcodes with N:", formatNumber(bpmap.readsWithN[i]) + " (" + getPercents(bpmap.readsWithN[i], bpmap.totalReads[i]) + "%)");
        for (int j=0; j<4; j++){        
            outputRow(ofs, polymsg + ATCG_BASES[j] + ":", formatNumber(bpmap.polyReads[i][j]) + " (" + getPercents(bpmap.polyReads[i][j], bpmap.totalReads[i]*bpmap.segment) + "%)");
        }
        outputRow(ofs, "duplicated barcodes:", formatNumber(bpmap.dupReads[i]) + " (" + getPercents(bpmap.dupReads[i], bpmap.totalReads[i]*bpmap.segment) + "%)");
        outputRow(ofs, "EST neighbor duplicated barcodes:", formatNumber(bpmap.ESTdupReads[i]) + " (" + getPercents(bpmap.ESTdupReads[i], bpmap.totalReads[i]*bpmap.segment) + "%)");
        outputRow(ofs, "Q10 bases:", getPercents(bpmap.readsQ10[i], bpmap.totalBase[i]) + "%");
        outputRow(ofs, "Q20 bases:", getPercents(bpmap.readsQ20[i], bpmap.totalBase[i]) + "%");
        outputRow(ofs, "Q30 bases:", getPercents(bpmap.readsQ30[i], bpmap.totalBase[i]) + "%");
        ofs << "</table>\n";
        ofs << "</div>\n";
    }
    /**
    if (mOptions->drawHeatMap.getQ10Tiff){
        string shrinkTiffName = basename(mOptions->drawHeatMap.q10TiffFile) + ".shrink.png";
        ofs << "<div class='subsection_title' onclick=showOrHide('q10Tiff')>chip q10 heatMap</div>\n";
        ofs << "<div id='q10Tiff'>\n";
        ofs << "<img src='" << shrinkTiffName << "' , alt='there is no image'>\n";
        ofs << "</div>\n";
        ofs << "<br>\n";
    }
    */
}

void HtmlReporter::printHeader(ofstream& ofs){
    ofs << "<html><head><meta http-equiv=\"content-type\" content=\"text/html;charset=utf-8\" />";
    ofs << "<title>spatialRNAextraction report " + getCurrentSystemTime() + " </title>";
    printJS(ofs);
    printCSS(ofs);
    ofs << "</head>";
    ofs << "<body><div id='container'>";
}

void HtmlReporter::printCSS(ofstream& ofs){
    ofs << "<style type=\"text/css\">" << endl;
    ofs << "td {border:1px solid #dddddd;padding:5px;font-size:12px;}" << endl;
    ofs << "table {border:1px solid #999999;padding:2x;border-collapse:collapse; width:800px}" << endl;
    ofs << ".col1 {width:240px; font-weight:bold;}" << endl;
    ofs << ".adapter_col {width:500px; font-size:10px;}" << endl;
    ofs << "img {padding:30px;}" << endl;
    ofs << "#menu {font-family:Consolas, 'Liberation Mono', Menlo, Courier, monospace;}" << endl;
    ofs << "#menu a {color:#0366d6; font-size:18px;font-weight:600;line-height:28px;text-decoration:none;font-family:-apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol'}" << endl;
    ofs << "a:visited {color: #999999}" << endl;
    ofs << ".alignleft {text-align:left;}" << endl;
    ofs << ".alignright {text-align:right;}" << endl;
    ofs << ".figure {width:800px;height:600px;}" << endl;
    ofs << ".header {color:#ffffff;padding:1px;height:20px;background:#000000;}" << endl;
    ofs << ".section_title {color:#ffffff;font-size:20px;padding:5px;text-align:left;background:#76b3df; margin-top:10px;}" << endl;
    ofs << ".subsection_title {font-size:16px;padding:5px;margin-top:10px;text-align:left;color:#76b3df}" << endl;
    ofs << "#container {text-align:center;padding:3px 3px 3px 10px;font-family:Arail,'Liberation Mono', Menlo, Courier, monospace;}" << endl;
    ofs << ".menu_item {text-align:left;padding-top:5px;font-size:18px;}" << endl;
    ofs << ".highlight {text-align:left;padding-top:30px;padding-bottom:30px;font-size:20px;line-height:35px;}" << endl;
    ofs << "#helper {text-align:left;border:1px dotted #fafafa;color:#777777;font-size:12px;}" << endl;
    ofs << "#footer {text-align:left;padding:15px;color:#ffffff;font-size:10px;background:#76b3df;font-family:Arail,'Liberation Mono', Menlo, Courier, monospace;}" << endl;
    ofs << ".kmer_table {text-align:center;font-size:8px;padding:2px;}" << endl;
    ofs << ".kmer_table td{text-align:center;font-size:8px;padding:0px;color:#ffffff}" << endl;
    ofs << ".sub_section_tips {color:#999999;font-size:10px;padding-left:5px;padding-bottom:3px;}" << endl;
    ofs << "</style>" << endl;
}

void HtmlReporter::printJS(ofstream& ofs){
    ofs << "<script src='https://cdn.plot.ly/plotly-latest.min.js'></script>" << endl;
    ofs << "\n<script type=\"text/javascript\">" << endl;
    ofs << "    function showOrHide(divname) {" << endl;
    ofs << "        div = document.getElementById(divname);" << endl;
    ofs << "        if(div.style.display == 'none')" << endl;
    ofs << "            div.style.display = 'block';" << endl;
    ofs << "        else" << endl;
    ofs << "            div.style.display = 'none';" << endl;
    ofs << "    }" << endl;
    ofs << "</script>" << endl;
}

const string HtmlReporter::getCurrentSystemTime()
{
  auto tt = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  struct tm* ptm = localtime(&tt);
  char date[60] = {0};
  sprintf(date, "%d-%02d-%02d      %02d:%02d:%02d",
    (int)ptm->tm_year + 1900,(int)ptm->tm_mon + 1,(int)ptm->tm_mday,
    (int)ptm->tm_hour,(int)ptm->tm_min,(int)ptm->tm_sec);
  return std::string(date);
}

void HtmlReporter::printFooter(ofstream& ofs){
    ofs << "\n<br>";
    ofs << "</div>" << endl;
    ofs << "<div id='footer'> ";
    //ofs << "<p>"<<command<<"</p>";
    ofs << "ST_BarcodeMap: getBarcode2Position map " << ", at " << getCurrentSystemTime() << " </div>";
    ofs << "</body></html>";
}

string HtmlReporter::list2string(long* list, int size) {
    stringstream ss;
    for (int i = 0; i < size; i++) {
        ss << list[i];
        if (i < size - 1)
            ss << ",";
    }
    return ss.str();
}
