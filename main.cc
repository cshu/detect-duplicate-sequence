#include <iostream>
#include <unordered_map>
#include <cassert>
#define NO_NEED_TO_AVERT_RACES
#include <cpprs.h>

using namespace std;

int main(int argc, char *argv[])
{
	try{
		if(argc<2)return 1;
		vector<vector<char>> fileco(argc-1);//optimize use make_unique<vector[]>(...)
		for(int i=argc-2;;--i){
			fileco[i]=readwholefileintovectorchar(argv[i+1]);
			fileco[i].shrink_to_fit();
			if(!i)break;
		}
#define siofblock 0x10//you need at least siofblock*2 consecutive bytes to be included in output, otherwise it's not fair (some 16-byte sequences are not searched, some 16-byte sequences are searched)
		unordered_map<string,vector<pair<int,size_t>>> sofblock;
		for(int i=argc-2;;--i){
			for(size_t doff=0;fileco[i].size()-doff>=siofblock;doff+=siofblock){
				auto pblock=fileco[i].data()+doff;
				for(const auto &elemofres:sofblock){
					//if(string::npos==elemofres.first.find(pblock,0,siofblock))goto nextblock;//note you can actually comment this line to include more extensive search
					auto lengthofblock=elemofres.first.size();
					for(const auto &occur:elemofres.second){//optimize no need to iterate all
						if(occur.first==i){
							if(occur.second<doff){
								if(doff-occur.second<lengthofblock)goto nextblock;
							}else{
								if(occur.second-doff<siofblock)goto nextblock;
							}
						}
					}
				}
				{
					vector<pair<int,size_t>> listofc;
					for(int j=argc-2;;--j){
						char *pch;
						size_t seoff=0;
						nextse:
						auto hsl=fileco[j].size()-seoff;
						if(hsl>=siofblock){
							UNSAFE_MEMMEMs(pch,fileco[j].data()+seoff,hsl,pblock,siofblock,{
								goto nextfile;
							},{
								seoff=pch-1-fileco[j].data();
								for(const auto &elemofres:sofblock){
									auto lengthofblock=elemofres.first.size();
									for(const auto &occur:elemofres.second){//optimize no need to iterate all
										if(occur.first==j){
											if(occur.second<seoff){
												if(seoff-occur.second<lengthofblock){++seoff;goto nextse;}
											}else{
												if(occur.second-seoff<siofblock){++seoff;goto nextse;}
											}
										}
									}
								}
								listofc.push_back(pair<int,size_t>(j,seoff));
								++seoff;
								goto nextse;
							})
						}
						nextfile:
						if(!j)break;
					}
					assert(listofc.size());
					if(listofc.size()==1){
						goto nextblock;
					}
					//all candidates are located, start expansion:
					//1. try to expand to left and right, as far as possible, ensuring all candidates are identical like    |<--[]-->|
					//2. if the new length is smaller than siofblock*2, use a sliding window like 0...|<--[]-->|...siofblock*2
					//   try all windows, and choose the window that ends up with highest frequency (without counting overlapped ones)
					//3. remove overlapped ones, and then try the second time to expand just like step 1
					//4. add to map
					size_t lenofcblock=siofblock;
#define expandbothsides(ui)\
					for(;;){\
						for(const auto &candidate:listofc)\
							if(!candidate.second)\
								goto ui##endofleftex;\
						auto testch=fileco[listofc[0].first][listofc[0].second-1];\
						for(const auto &candidate:listofc)/*//optimize skip first iteration*/\
							if(testch!=fileco[candidate.first][candidate.second-1])\
								goto ui##endofleftex;\
						for(const auto &candidate:listofc){\
							for(const auto &elemofres:sofblock){\
								auto lengthofblock=elemofres.first.size();\
								for(const auto &occur:elemofres.second){/*//optimize skip first iteration*/\
									if(occur.first==candidate.first)\
										if(occur.second+lengthofblock==candidate.second)\
											goto ui##endofleftex;\
								}\
							}\
						}\
						for(auto &candidate:listofc)\
							--candidate.second;\
						++lenofcblock;\
					}\
					ui##endofleftex:\
					for(;;){\
						for(const auto &candidate:listofc)\
							if(lenofcblock==fileco[candidate.first].size()-candidate.second)\
								goto ui##endofrightex;\
						auto testch=fileco[listofc[0].first][listofc[0].second+lenofcblock];\
						for(const auto &candidate:listofc)/*//optimize skip first iteration*/\
							if(testch!=fileco[candidate.first][candidate.second+lenofcblock])\
								goto ui##endofrightex;\
						for(const auto &candidate:listofc){\
							for(const auto &elemofres:sofblock){\
								for(const auto &occur:elemofres.second){/*//optimize skip first iteration*/\
									if(occur.first==candidate.first)\
										if(occur.second==candidate.second+lenofcblock)\
											goto ui##endofrightex;\
								}\
							}\
						}\
						++lenofcblock;\
					}\
					ui##endofrightex:;
expandbothsides(first_ex_)
					if(lenofcblock<siofblock*2){
						vector<pair<int,size_t>> bestwindow;
						for(auto lenofrecession=siofblock*2-lenofcblock;;){//slide window
							vector<decltype(listofc)> currentwindow{listofc};
							for(auto iter=currentwindow[0].cbegin();iter!=currentwindow[0].cend();){
								if(lenofrecession>iter->second || siofblock*2-lenofcblock-lenofrecession>fileco[iter->first].size()-iter->second)
									iter=currentwindow[0].erase(iter);
								else
									++iter;
							}
							for(auto &candidate:currentwindow[0])
								candidate.second-=lenofrecession;
							for(auto iter=currentwindow[0].cbegin();iter!=currentwindow[0].cend();){
								for(const auto &elemofres:sofblock){
									auto lengthofblock=elemofres.first.size();
									for(const auto &occur:elemofres.second){//optimize no need to iterate all
										if(occur.first==iter->first)
											if(check2intervalsoverlap(occur.second,occur.second+lengthofblock,iter->second,iter->second+siofblock*2))goto erasei;
									}
								}
								++iter;
								continue;
								erasei:
								iter=currentwindow[0].erase(iter);
							}
							for(auto iter=currentwindow[0].cbegin();iter!=currentwindow[0].cend();){
								if(memcmp(fileco[currentwindow[0][0].first].data()+currentwindow[0][0].second,fileco[iter->first].data()+iter->second,siofblock*2)){
									for(auto &wind:currentwindow){
										if(!memcmp(fileco[wind[0].first].data()+wind[0].second,fileco[iter->first].data()+iter->second,siofblock*2)){
											wind.push_back(*iter);
											goto eraserelocated;
										}
									}
									currentwindow.push_back(vector<pair<int,size_t>>{*iter});
									eraserelocated:
									iter=currentwindow[0].erase(iter);
								}else{
									++iter;
								}
							}
							for(auto &wind:currentwindow){
								int previousfirst=-1;
								size_t previoussecond;
								for(auto windbl=wind.cbegin();windbl!=wind.cend();){
									if(windbl->first==previousfirst){
										assert(windbl->second>previoussecond);
										if(siofblock*2>windbl->second-previoussecond){
											windbl=wind.erase(windbl);
											continue;
										}
									}
									previousfirst=windbl->first;
									previoussecond=windbl->second;
									++windbl;
								}
							}
							for(auto &wind:currentwindow){//todo when you have multple groups, retain all of them rather than retain only the biggest one?
								if(wind.size()>bestwindow.size())
									bestwindow=move(wind);
							}
							--lenofrecession;
							if(!lenofrecession)break;
						}
						if(bestwindow.size()<2)goto nextblock;
						listofc=move(bestwindow);
						lenofcblock=siofblock*2;
expandbothsides(second_ex_)
					}
					//undone trim whitespace characters, and then check length once more. if too short, skip
					sofblock.emplace(string(fileco[listofc[0].first].data()+listofc[0].second,lenofcblock),move(listofc));
				}
				nextblock:;
			}
			if(!i)break;
		}
		for(const auto &block:sofblock){
			cout<<"SEQUENCE BEGIN\n"<<block.first<<"\nSEQUENCE END\n";
			for(const auto &occur:block.second){
				cout<<argv[occur.first+1]<<' '<<occur.second<<'\n';
			}
		}
		return 0;
	}catch(const std::exception &e){
		STD_CLOG_FILE_FUNC_LINE_EX_FLUSH_NOEXCEPTs(e)
	}catch(...){
		STD_CLOG_FILE_FUNC_LINE_FLUSH_NOEXCEPTs
	}
	return 1;
}
