% Cell interaction
 % In this work we consider the mutal  interactions between the cells
 % like what 
 % we have in the bilogical systems. In this model the release rate fuction
 % comrise a part from the natural realease and a part as a feedback from the
 % target cell.
clc;
clear all
close all;

%load the simulation parameters
Parameters

% Spherical bessel function
SphBess=@(x,n)((pi/2./x).^.5.*besselj(n+0.5,x));
SphHank=@(x,n)((pi/2./x).^.5.*besselh(n+0.5,2,x));
SphHankd=@(x,n)(2*n+1).^-1.*(n.*SphHank(x,n-1)-(n+1).*SphHank(x,n+1));
SphBessd=@(x,n)(2*n+1).^-1.*(n.*SphBess(x,n-1)-(n+1).*SphBess(x,n+1));
SphBessn=@(x,n)((pi/2./x).^.5.*bessely(n+0.5,x));
SphBessdn=@(x,n)(2*n+1).^-1.*(n.*SphBessn(x,n-1)-(n+1).*SphBessn(x,n+1));



RR_I=zeros(length(d),2*length(t0)+1); %number of released vesicles by the immune cell

j=1;
% Release coefficient of tumor EV


i=1;
% distanse between cell membranes
for rm=d %[output:group:68d6f3f8]

    %counter
    (rm-d(1))/(d(end)-d(1))*100 %[output:67255e37] %[output:78d58831] %[output:7947801a] %[output:3f93e6ba] %[output:1fecd566] %[output:1927dbe6] %[output:43f818e8] %[output:3cedce63] %[output:0512aa81] %[output:8eeb2360] %[output:9da3be31] %[output:21d48ef2] %[output:4909baa8] %[output:5174914a] %[output:6850a378] %[output:062a608f] %[output:93e9ac50] %[output:2bc75bcf] %[output:3e15590b] %[output:5ad63e60] %[output:2f96b81a] %[output:4ec6be17] %[output:4329c6a7] %[output:1e6b9578] %[output:8c461895] %[output:4564161c] %[output:399942c0] %[output:52568e74] %[output:9c1b4581] %[output:23eb66a0] %[output:88845499] %[output:56e4ebf8] %[output:674c21b1] %[output:50180601] %[output:4425a3f5] %[output:9e0f0288] %[output:8cfeca85] %[output:72be3094] %[output:9f7ff3e5] %[output:65aff999] %[output:66ecea3f] %[output:8b3e0ea9] %[output:1e9f7840] %[output:29c58d41] %[output:799d7fac] %[output:70b47185] %[output:0a877eae] %[output:8623ed0d] %[output:71db15d8] %[output:97bd5390] %[output:4c04bdb9] %[output:3fd60ce3] %[output:443f10fc] %[output:5937b00c] %[output:890d5137] %[output:8359d3ef] %[output:8f02e387] %[output:21792d8c] %[output:139229e6] %[output:9093d9a2] %[output:0bf9fdad] %[output:3c14fac7] %[output:43fb5f33] %[output:1b4ae426] %[output:2c9c9ab7] %[output:6a73c26d] %[output:7102f20f] %[output:88c830a0] %[output:905ae766] %[output:0b267d00] %[output:1574076d] %[output:5d22b67e] %[output:76bd1a6c] %[output:9d410eca] %[output:384efeb9] %[output:1252cc27] %[output:76bb0a29] %[output:1b65c340] %[output:328d6f15] %[output:52181767] %[output:9f724b22] %[output:4cb6043d] %[output:33c84327] %[output:1568c1f2] %[output:2a92733c] %[output:45aa9910] %[output:9e61474c] %[output:36bb927b] %[output:8999f1a0] %[output:64e8eaab] %[output:50fee51b] %[output:063188c6] %[output:1f9feb4f] %[output:9ce6d470] %[output:599f3320] %[output:8b89ab8d] %[output:011a88b0] %[output:3214ef59] %[output:6cc52d94] %[output:68743e90] %[output:78c8c4cd] %[output:07d9edb2] %[output:22cd13c0] %[output:43fc3e51] %[output:9c915a32] %[output:96359e81] %[output:38d1403a] %[output:33658f73] %[output:0e7778ec] %[output:0d19b4f5] %[output:9982c429] %[output:689687c7] %[output:9a4a49da] %[output:7cfd1787] %[output:5e63efe4] %[output:223a66be] %[output:6a935bea] %[output:3371492d] %[output:870790e4] %[output:34538248] %[output:3df6f7f0] %[output:28016d60] %[output:1aadc492] %[output:694a84f7] %[output:13ce17c3] %[output:37358051] %[output:7af49e47] %[output:5bda90ad] %[output:02d3e201] %[output:6f1db92d] %[output:91134c5b] %[output:59afc586] %[output:22f9fa0e] %[output:9da7cd9c] %[output:1b83bb7e] %[output:4751edfe] %[output:024c3c35] %[output:4c3eb3b0] %[output:51756178] %[output:4b09fc72] %[output:8dba2d56] %[output:139759b8] %[output:8ccb9a5c] %[output:99724e8f] %[output:1ae99635] %[output:60393c8d] %[output:7d5728bd] %[output:54d1675b] %[output:069a3da4] %[output:83ff8d41] %[output:5d696025] %[output:01b69dd4] %[output:535f2d3e] %[output:75e40f95] %[output:78200454] %[output:76d7aefc] %[output:215e545d] %[output:2e22e3f5] %[output:03de8895] %[output:09cba00b] %[output:7b65ba6e] %[output:77b435d0] %[output:75272ce8] %[output:84a92e87] %[output:2d215bbb] %[output:0c4442d6] %[output:964605ac] %[output:6aacb034] %[output:7a8c0d75] %[output:6555e12f] %[output:6c1e140e] %[output:9c9cc535] %[output:7bd66a30] %[output:5b9e7555] %[output:1d98a0ea] %[output:5fd96aae] %[output:706ec828] %[output:9b4de024] %[output:2c02632f] %[output:226614ba] %[output:514b8679] %[output:5a0921bd] %[output:6eb7f962] %[output:25f9a8f7] %[output:1fcd2943] %[output:613dfbc8] %[output:6fb2c1d3] %[output:6a77544f] %[output:48d54184] %[output:1ab92c76] %[output:2b17c194] %[output:3d6322a2] %[output:4eba9adf] %[output:01fea283] %[output:84bda8bb] %[output:88b358bf] %[output:80fd3a13] %[output:9a8f6e32] %[output:0a4efcb4] %[output:81bec6cd] %[output:18be26ce] %[output:9656efa4] %[output:1b67f9a1] %[output:2bfe19c1] %[output:8a3f7eb6] %[output:38b7a692] %[output:2929b9e0] %[output:1880a2b7] %[output:8211b003] %[output:1021cb98] %[output:136342e5] %[output:00c8a90d] %[output:5bc653c5] %[output:52e165b8] %[output:69718805] %[output:4454b93e] %[output:2eb0fc25] %[output:2f1ebc73] %[output:6ac1855c] %[output:46048279] %[output:5a9134e4] %[output:75bec839] %[output:266dd468] %[output:7b45b5e8] %[output:4c6b05e4] %[output:6a47fc37] %[output:8c1ecd7e] %[output:363e31cc] %[output:165731f1] %[output:0eeb0472] %[output:1b195a90] %[output:270a5f25] %[output:44414d94] %[output:23330413] %[output:3a6cb781] %[output:95fbffbe] %[output:64923a30] %[output:99673cf8] %[output:0c699fe3] %[output:897f6040] %[output:9cdbd7e5] %[output:35624b29] %[output:23ecd33e] %[output:030acaa6] %[output:54a60db0] %[output:95e948af] %[output:45fa1c6c] %[output:58cff72d] %[output:274d8d71] %[output:66668248] %[output:07fe4094] %[output:5bbb967a] %[output:4d25b75b] %[output:603fdb6a] %[output:955a827a] %[output:755e79a7] %[output:40e509b8] %[output:321f1155] %[output:6f92ed21] %[output:7b672439] %[output:1ab9803a] %[output:149db7da] %[output:0b5b4b5a] %[output:1c31cd71] %[output:2f6f7c77] %[output:22091076] %[output:47d5fa64] %[output:74c79afc] %[output:96ff3f78] %[output:801a9c5b] %[output:132d41c4] %[output:81cc19d7] %[output:1696d515] %[output:37fb1110] %[output:31aa9287] %[output:937994ee] %[output:7b9b84f7] %[output:78338b7b] %[output:8ea3a2c2] %[output:8c20c088] %[output:7fb82433] %[output:396ebf35] %[output:7e334680] %[output:04c5d578] %[output:50b9dcd4] %[output:0c3c5267] %[output:82408e79] %[output:381f08aa] %[output:65403290] %[output:059555df] %[output:5b3c959b] %[output:110c4a7c] %[output:6fcf1ccf] %[output:545316d8] %[output:4e8e2520] %[output:44fb241d] %[output:71148cef] %[output:3b4a7852] %[output:6823c738] %[output:7a0d0ca7] %[output:9ef39bf5] %[output:51d8c4b6] %[output:6192d30b] %[output:3c7ccce6] %[output:177d8c00] %[output:356bcfe5] %[output:28080b05] %[output:9a14a5bb] %[output:28ade82f] %[output:1a641902] %[output:28593b5b] %[output:70c22289] %[output:967637ea] %[output:7a96a34e] %[output:7885d04c] %[output:17e623ac] %[output:9d4942fb] %[output:46c0b5dd] %[output:39dccf90] %[output:8fe6fe1b] %[output:2eda21b3] %[output:71672212] %[output:69273952] %[output:3f6f2ce7] %[output:058742a4] %[output:431c0cb8] %[output:53b6ab01] %[output:27cafebb] %[output:18da230d] %[output:6fea10df] %[output:96d15cf4] %[output:53c24019] %[output:7f702b9d] %[output:28180a47] %[output:35aa60ec] %[output:9661b423] %[output:9f1422bf] %[output:6754ab06] %[output:956eeda1] %[output:7651eee1] %[output:92380a8b] %[output:6be12667] %[output:60599c40] %[output:78e0c533] %[output:84ddf029] %[output:4b178b3d] %[output:8863680a] %[output:1d45b848] %[output:3639c26b] %[output:87831e8d] %[output:8477d437] %[output:4a3082f5] %[output:95e26c08] %[output:7004646e] %[output:11612855] %[output:71305314] %[output:06f8ded7] %[output:1a641fef] %[output:569539f6] %[output:264acde4] %[output:76b6f09b] %[output:620bc512] %[output:7da73d49] %[output:1885aa20] %[output:0c2fc8e8] %[output:7d4eaf7c] %[output:6563ea6b] %[output:2a17c920] %[output:0314a60f] %[output:6a43efc5] %[output:2431247b] %[output:1b85c110] %[output:744a84f8] %[output:6205f959] %[output:8407e28e] %[output:9ea3a873] %[output:24765adb] %[output:019b01a0] %[output:5bf46f37] %[output:58acdc61] %[output:78e17c4e] %[output:26c3f164] %[output:8257d51a] %[output:41e6b9b2] %[output:93f978bc] %[output:99db324d] %[output:838e2deb] %[output:489ee308] %[output:381d25ba] %[output:3732681b] %[output:395eedaa] %[output:9f2efd6d] %[output:0b064d9d] %[output:0c3b85e8] %[output:89b380fd] %[output:0fca8b88] %[output:228e289e] %[output:6943a535] %[output:305a5172] %[output:1ee3b48e] %[output:56b344e8] %[output:5f4d8a5b] %[output:964b5714] %[output:9dbbea33] %[output:10a0b4f5] %[output:760d7c0c] %[output:21f6e935] %[output:69df8682] %[output:79e75d23] %[output:5b6d7efd] %[output:85a938c2] %[output:1ef591da] %[output:324f5a68] %[output:05d2a5cf] %[output:82d6592f] %[output:255ad6d8] %[output:5289743b] %[output:85401496] %[output:40b738d3] %[output:1bd329c7] %[output:3763fbd1] %[output:96ee24d6] %[output:99e01ee0] %[output:3bc25895] %[output:1baae318] %[output:0c76f1d0] %[output:4945bf43] %[output:96e18e9a] %[output:14ec4758] %[output:69f6ec07] %[output:9915aa92] %[output:2060f380] %[output:5e2902ee] %[output:1406a359] %[output:564aba1d] %[output:3e3baa49] %[output:548e9430] %[output:7ad947db] %[output:3b460793] %[output:84f61f1f] %[output:151f66f7] %[output:9af9a3a3] %[output:191448a2] %[output:9db83ffc] %[output:7127c4b6] %[output:00306098] %[output:0b054eb4] %[output:0e436b1d] %[output:01ba0877] %[output:8ed2f7b4] %[output:072b8c96] %[output:417e144e] %[output:9085762e] %[output:483f17bc] %[output:34efc0bc] %[output:47c8f7d7] %[output:9d56ccc0] %[output:9035a2ce] %[output:3f05b86c] %[output:3f61d08d] %[output:10961157] %[output:5914a913] %[output:5118c604] %[output:48cbe5d4] %[output:26829e64] %[output:27b587b6] %[output:380980d0] %[output:6e262f29] %[output:6defad62] %[output:4cec0fec] %[output:0b28b513] %[output:407931e4] %[output:0cff5c09] %[output:2d32fbf1] %[output:33b53d72] %[output:6e4430e8] %[output:2fac5439] %[output:463cf34f] %[output:77639cee] %[output:736f75a4] %[output:37c8c8ff] %[output:8fcdcf3d] %[output:0cbb6239] %[output:70da4a1f] %[output:8c72b683] %[output:7732bf28] %[output:38309d76] %[output:121ab693] %[output:72479b01] %[output:144114a2] %[output:89768699] %[output:859e995f] %[output:5340c61d] %[output:259b9533] %[output:4d582263] %[output:12a6e5ac] %[output:5b33f8ce] %[output:804edcbe] %[output:04cd8b18] %[output:8ac5e673] %[output:7f2f9789] %[output:78836302] %[output:1d244cce] %[output:3cae733e] %[output:5a98d013] %[output:4a2f020b] %[output:3ac4683b] %[output:5e401a06] %[output:1cc3633b] %[output:970ad9af] %[output:1950352c] %[output:3cddb129] %[output:120fe386] %[output:05c29855] %[output:56c1b452] %[output:20072cf5] %[output:45194e1b] %[output:19b62263] %[output:994f6dfe] %[output:426eea91] %[output:47eb7e26] %[output:4dba9289] %[output:7b54bd6e] %[output:32fa3add] %[output:7046e0a8] %[output:8dea4457] %[output:18c8dcef] %[output:9da87ee6] %[output:6562bce7] %[output:729b5839] %[output:51d0d79b] %[output:84a8bada] %[output:94ced4b9] %[output:4a908277] %[output:533db7d4] %[output:7197219b] %[output:0314c648] %[output:7e27bbaa] %[output:3843cfd1] %[output:01e79ea5] %[output:0ae75e26] %[output:5ceb97da] %[output:2852c5ef] %[output:9511a384] %[output:5cbb6b72] %[output:20bb45ee] %[output:79b6ff1a] %[output:36dc67ec] %[output:691d1057] %[output:36de0900] %[output:6098be5d] %[output:6759a4b1] %[output:4056dafc] %[output:056303c7] %[output:49cd1bbb] %[output:9995b870] %[output:0e356ab9] %[output:8f00889b] %[output:74e19455] %[output:8ec79a8c] %[output:28fea522] %[output:365f4b1e] %[output:3d47d553] %[output:98d7f80d] %[output:5ab8911d] %[output:26d022a8] %[output:6e63325a] %[output:9bd48408] %[output:7d06c9f5] %[output:63df2a33] %[output:6d414100] %[output:1c072ab4] %[output:229ec812] %[output:41a392e7] %[output:178c0d9d] %[output:06d54718] %[output:2f6077ad] %[output:182a79a5] %[output:22076f86] %[output:5b538180] %[output:78a90a2a] %[output:8e92e969] %[output:1d3e6eb6] %[output:3e215f15] %[output:2fbbf178] %[output:0fc2bba4] %[output:025a2069] %[output:95010784] %[output:7467d443] %[output:25231ef8] %[output:1101d565] %[output:778212f0] %[output:5e610b10] %[output:0722b6de] %[output:4861de90] %[output:14f02c7c] %[output:7d88f410] %[output:586ac50b] %[output:6f64eba1] %[output:80223990] %[output:13712d29] %[output:0ec1931e] %[output:81dd539c] %[output:9fac415c] %[output:5f116e0e] %[output:1dc8fc8b] %[output:23e2e1b0] %[output:1f609d46] %[output:673028df] %[output:9d2c4fa0] %[output:178db462] %[output:9a6c55d8] %[output:2e105267] %[output:8cff8852] %[output:421cc68b] %[output:2a0d85a2] %[output:3353f2fd] %[output:664c92a7] %[output:7e658694] %[output:15162885] %[output:51ad56c2] %[output:55c1b230] %[output:81da06be] %[output:957fabd9] %[output:3cd8fc4b] %[output:9e23b446] %[output:361f5440] %[output:2215f5c9] %[output:2ff88c13] %[output:293ebc1d] %[output:5438b0db] %[output:8b2ee12b] %[output:9714b09f] %[output:424feca3] %[output:21834e3e] %[output:5b03f33d] %[output:341272de] %[output:01292d95] %[output:8486e9a9] %[output:51c99304] %[output:9e2d8635] %[output:2657df05] %[output:2f8562f6] %[output:239b917a] %[output:3830e1a1] %[output:2d8d1fa7] %[output:0693e091] %[output:964cc592] %[output:752a08cc] %[output:8c6c6a54] %[output:614c5ccf] %[output:08907412] %[output:6a20dbe7] %[output:68113b07] %[output:2eef1172] %[output:56d33e48] %[output:69be270b] %[output:5d3a4ddc] %[output:83a85315] %[output:7e321cb0] %[output:0995325f] %[output:9df4b145] %[output:2b494034] %[output:4c5520cd] %[output:96ca4210] %[output:2c03b7e6] %[output:8d0b49d5] %[output:7d8b3732] %[output:43159cbb] %[output:7e55382e] %[output:12eb3bd4] %[output:66b71042] %[output:05c7e8e9] %[output:98c9f9f4] %[output:8fc83dc5] %[output:62923a04] %[output:07e40f0f] %[output:73cc72ce] %[output:16ae563c] %[output:7ece0681] %[output:21960f4f] %[output:8e0387f5] %[output:464479c2] %[output:047489c9] %[output:94f83b2c] %[output:181206d1] %[output:0100eb18] %[output:05ab07ae] %[output:036d51ec] %[output:5a4d8a80] %[output:7511b3b3] %[output:3820b7f8] %[output:29753cd6] %[output:80166951] %[output:6de4c5de] %[output:2db7f1f8] %[output:68c14a0b] %[output:3fdf370d] %[output:7e55fb39] %[output:8c865a76] %[output:15fe009f] %[output:874066da] %[output:94230be8] %[output:8152fc03] %[output:4ab10d19] %[output:4d24f2f4] %[output:20699bd6] %[output:1fac4615] %[output:078bc257] %[output:087401ef] %[output:90aebbfb] %[output:8e9deb37] %[output:3d474907] %[output:2f4a9930] %[output:6ca79501] %[output:7702b830] %[output:97a66719] %[output:0d0b52fd] %[output:228b9af6] %[output:14bb6bed] %[output:1e94472b] %[output:257b54dd] %[output:0f3283fb] %[output:870188e8] %[output:393d000d] %[output:6a0f3061] %[output:702b1eef] %[output:0fb12f55] %[output:7748a852] %[output:19f568c0] %[output:3a60f553] %[output:553f48f4] %[output:2b4d7004] %[output:45a2957e] %[output:621a2a92] %[output:42454fd2] %[output:68e11212] %[output:5692ca9b] %[output:85e293fc] %[output:72da5d47] %[output:746fbf2e] %[output:425da69e] %[output:44467add] %[output:2eb9db58] %[output:638663c1] %[output:959dce23] %[output:98df5fbd] %[output:2facbc44] %[output:447d832e] %[output:5f4bb0cb] %[output:0dbdb5cb] %[output:3d0925f0] %[output:48702547] %[output:3d0f8974] %[output:8ca5d660] %[output:1e9e16c4] %[output:7f057364] %[output:5f66681d] %[output:0b5835b3] %[output:30f3acd8] %[output:6e29d86c] %[output:6ef18a41] %[output:322bf010] %[output:6c78e3aa] %[output:68defd27] %[output:98b53265] %[output:09921628] %[output:9d0ca7d6] %[output:607e7357] %[output:0f9687b7] %[output:3afab700] %[output:1d0a6ebf] %[output:8a259384] %[output:507b6d0e] %[output:35f56edf] %[output:6f7ece41] %[output:0ba0aa8e] %[output:382e5f72] %[output:669b5316] %[output:52557c3b] %[output:0b941c1e] %[output:1def1396] %[output:0d5ec1f4] %[output:37f2f837] %[output:50a9f906] %[output:1eddf2b9] %[output:4647a7a7] %[output:4ff59d3b] %[output:655da3ac] %[output:478093d0] %[output:4151fff2] %[output:365b3483] %[output:43a15dcf] %[output:0b62e4d8] %[output:22b01d39] %[output:2ecf2c42] %[output:2925a644]

    % distanse between cell centers
    rc=rm+a_T+a_I;
    r_T=a_I+rm;
    r_I=a_T+rm;
    %val = @(x) a1*exp(-((x/3600-b1)/d1).^2) + a2*exp(-((x/3600-b2)/d2).^2);

    %% Evaluation of the CIR at the immune and tumor cells
    jj=1;
    tic

    % Frequency response for the tumor to immune cell H_TI(Natural release)
    for omega=pi*omega0

        Dw=D_T*(1i*omega+kr_I+ki_I);
        Lw=(1i*omega+ki_I)*kb_I;
        kdp=-(kd_T+1i*omega)/D_T;
        D=D_T;
        % Natural release
        r0=rc;
        a=a_I;
        l=0;
        A1=Dw.*sqrt(kdp)*SphBessdn(sqrt(kdp)*a,l)-Lw.*SphBessn(sqrt(kdp)*a,l);
        B1=Dw*sqrt(kdp)*SphBessd(sqrt(kdp)*a,l)-Lw.*SphBess(sqrt(kdp)*a,l);
        A2=SphBessn(sqrt(kdp)*r0,l);
        B2=SphBess(sqrt(kdp)*r0,l);
        C2=-SphHank(sqrt(kdp)*r0,l);
        A3=r0.^2*sqrt(kdp)*SphBessdn(sqrt(kdp)*r0,l);
        B3=r0.^2*sqrt(kdp)*SphBessd(sqrt(kdp)*r0,l);
        C3=-r0.^2*sqrt(kdp)*SphHankd(sqrt(kdp)*r0,l);
        B0=1/D./(-A3.*B1./A1+B3-C3.*(-A2.*B1./A1+B2)./C2);
        A0=-B1./A1.*B0;
        CR=(A0.*SphBessn(sqrt(kdp)*a,0)+B0.*SphBess(sqrt(kdp)*a,0))*1;

        N_bind(jj)=kb_I*CR./(1j*omega+kr_I+ki_I)*a^2;
        if isnan(N_bind(jj))
            N_bind(jj) =  N_bind(jj-1);
        end
        % Internalization frequency response for natural release
        H_TIn(jj)=ki_I*N_bind(jj);
        jj=jj+1;

    end
    jj=1;
    % Frequency response for the tumor to immune cell H_TI(Mutual release)
    for omega=pi*omega0

        Dw=D_T*(1i*omega+kr_I+ki_I);
        Lw=(1i*omega+ki_I)*kb_I;
        kdp=-(kd_T+1i*omega)/D_T;
        D=D_T;
        % Natural release
        r0=r_T;
        a=a_I;
        l=0;
        A1=Dw.*sqrt(kdp)*SphBessdn(sqrt(kdp)*a,l)-Lw.*SphBessn(sqrt(kdp)*a,l);
        B1=Dw*sqrt(kdp)*SphBessd(sqrt(kdp)*a,l)-Lw.*SphBess(sqrt(kdp)*a,l);
        A2=SphBessn(sqrt(kdp)*r0,l);
        B2=SphBess(sqrt(kdp)*r0,l);
        C2=-SphHank(sqrt(kdp)*r0,l);
        A3=r0.^2*sqrt(kdp)*SphBessdn(sqrt(kdp)*r0,l);
        B3=r0.^2*sqrt(kdp)*SphBessd(sqrt(kdp)*r0,l);
        C3=-r0.^2*sqrt(kdp)*SphHankd(sqrt(kdp)*r0,l);
        B0=1/D./(-A3.*B1./A1+B3-C3.*(-A2.*B1./A1+B2)./C2);
        A0=-B1./A1.*B0;
        CR=(A0.*SphBessn(sqrt(kdp)*a,0)+B0.*SphBess(sqrt(kdp)*a,0))*1;
        N_bind(jj)=kb_I*CR./(1j*omega+kr_I+ki_I)*a^2;
        if isnan(N_bind(jj))
            N_bind(jj) = N_bind(jj-1);
        end
        % Internalization frequency response for mutual release
        H_TIm(jj)=ki_I* N_bind(jj);
        jj=jj+1;

    end

    jj=1;

    % Frequency response for the immune to tumor cell H_IT(Natural release)
    for omega=pi*omega0

        Dw=D_I*(1i*omega+kr_T+ki_T);
        Lw=(1i*omega+ki_T)*kb_T;
        kdp=-(kd_I+1i*omega)/D_I;
        D=D_I;
        % Natural release
        r0=rc;
        a=a_T;
        l=0;
        A1=Dw.*sqrt(kdp)*SphBessdn(sqrt(kdp)*a,l)-Lw.*SphBessn(sqrt(kdp)*a,l);
        B1=Dw*sqrt(kdp)*SphBessd(sqrt(kdp)*a,l)-Lw.*SphBess(sqrt(kdp)*a,l);
        A2=SphBessn(sqrt(kdp)*r0,l);
        B2=SphBess(sqrt(kdp)*r0,l);
        C2=-SphHank(sqrt(kdp)*r0,l);
        A3=r0.^2*sqrt(kdp)*SphBessdn(sqrt(kdp)*r0,l);
        B3=r0.^2*sqrt(kdp)*SphBessd(sqrt(kdp)*r0,l);
        C3=-r0.^2*sqrt(kdp)*SphHankd(sqrt(kdp)*r0,l);
        B0=1/D./(-A3.*B1./A1+B3-C3.*(-A2.*B1./A1+B2)./C2);
        A0=-B1./A1.*B0;
        CR=(A0.*SphBessn(sqrt(kdp)*a,0)+B0.*SphBess(sqrt(kdp)*a,0))*1;
        N_bind(jj)=kb_T*CR./(1j*omega+kr_T+ki_T)*a^2;
        if isnan(N_bind(jj))
            N_bind(jj) = N_bind(jj-1);
        end
        % Internalization frequency response for natural release
        H_ITn(jj)=ki_T* N_bind(jj);
        jj=jj+1;

    end
    jj=1;
    % Frequency response for the tumor to immune cell H_TI(Mutual release)
    for omega=pi*omega0
        Dw=D_I*(1i*omega+kr_T+ki_T);
        Lw=(1i*omega+ki_T)*kb_T;
        kdp=-(kd_I+1i*omega)/D_I;
        D=D_I;
        % Natural release
        r0=r_I;
        a=a_T;
        l=0;
        A1=Dw.*sqrt(kdp)*SphBessdn(sqrt(kdp)*a,l)-Lw.*SphBessn(sqrt(kdp)*a,l);
        B1=Dw*sqrt(kdp)*SphBessd(sqrt(kdp)*a,l)-Lw.*SphBess(sqrt(kdp)*a,l);
        A2=SphBessn(sqrt(kdp)*r0,l);
        B2=SphBess(sqrt(kdp)*r0,l);
        C2=-SphHank(sqrt(kdp)*r0,l);
        A3=r0.^2*sqrt(kdp)*SphBessdn(sqrt(kdp)*r0,l);
        B3=r0.^2*sqrt(kdp)*SphBessd(sqrt(kdp)*r0,l);
        C3=-r0.^2*sqrt(kdp)*SphHankd(sqrt(kdp)*r0,l);
        B0=1/D./(-A3.*B1./A1+B3-C3.*(-A2.*B1./A1+B2)./C2);
        A0=-B1./A1.*B0;
        CR=(A0.*SphBessn(sqrt(kdp)*a,0)+B0.*SphBess(sqrt(kdp)*a,0))*1;
        N_bind(jj)=kb_T*CR./(1j*omega+kr_T+ki_T)*a^2;
        if isnan(N_bind(jj))
            N_bind(jj) = N_bind(jj-1);
        end
        % Internalization frequency response for mutual release
        H_ITm(jj)=ki_T* N_bind(jj);
        jj=jj+1;

    end

    % Closed loop Tumor rate in frequency
    % R_Tf=[0 R_Tf*max(omega0) conj(fliplr(R_Tf*max(omega0)))];
    % R_If=[0 R_If*max(omega0) conj(fliplr(R_If*max(omega0)))];
    % Omg=[domega*pi omega0*pi omega0*pi+pi*omega0(end)];

    H_ITn1=[0 H_ITn conj(fliplr(H_ITn))];
    H_ITm1=[0 H_ITm conj(fliplr(H_ITm))];
    H_TIn1=[0 H_TIn conj(fliplr(H_TIn))];
    H_TIm1=[0 H_TIm conj(fliplr(H_TIm))];

    Omg=[ 0 omega0*pi max(omega0*pi)+omega0*pi];

    % RR_Tf=(R_Tf+kre_T*H_ITn.*R_If)./(1-kre_T*kre_I.*H_ITm.*H_TIm);
    % Internalized EVs rate in Tumor due to natural release from Immune
    I_nTf=H_ITn1.*R_If;
    % Internalized EVs rate in Immune due to natural release from Tumor
    I_nIf=H_TIn1.*R_Tf;
    % Overall internalization rate  fuction at Tumor
    I_Tf=(H_ITn1.*R_If+kre_I*H_ITm1.*H_TIn1.*R_Tf)./...
        (1-kre_T*kre_I.*H_TIm1.*H_ITm1);
    % I_Tf=(1j*Omg.*H_ITn.*R_If+kre_I*H_ITm.*H_TIn.*R_Tf)./...
    %     (-Omg.^2-kre_T*kre_I.*H_TIm.*H_ITm);
    % Overall internalization rate  fuction at Immune
    I_If=(H_TIn1.*R_Tf+kre_T*H_TIm1.*H_ITn1.*R_If)./...
        (1-kre_I*kre_T.*H_ITm1.*H_TIm1);
    % Overall tumor release
    RR_Tf=R_Tf+kre_T.*I_Tf;
    % Overall immune release
    RR_If=R_If+kre_I.*I_If;

    sys=frd(I_If,Omg);

    % Tumor rate in time the domain
    RR_T=real(ifft(RR_Tf))*max(omega0);
    RR_I(i,:)=real(ifft(RR_If))*max(omega0);
    R_T=real(ifft(R_Tf))*max(omega0);
    H_ITnt=real(ifft(H_ITn))*max(omega0);
    H_TInt=real(ifft(H_TIn))*max(omega0);
    H_ITmt=real(ifft(H_ITm))*max(omega0);
    H_TImt=real(ifft(H_TIm))*max(omega0);
    I_T=real(ifft(I_Tf));
    I_I=real(ifft(I_If));
    I_nT=real(ifft(I_nTf))*max(omega0);
    IR(i,j)=max(RR_I(i,:));
    i=i+1;
end %[output:group:68d6f3f8]
%%
%[text] #### Saving the Dataset
%the first column stores the predictor (distance between the cells)
%the next column stores the amount of vesicles with time along each row

Dataset=array2table(horzcat(d',(RR_I(:,1:floor(end/2)))));

save Dataset_cell2cell dist slope peak

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":38.9}
%---
%[output:67255e37]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"0"},"version":0}
%---
%[output:78d58831]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"0.1250"},"version":0}
%---
%[output:7947801a]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"0.2500"},"version":0}
%---
%[output:3f93e6ba]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"0.3750"},"version":0}
%---
%[output:1fecd566]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"0.5000"},"version":0}
%---
%[output:1927dbe6]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"0.6250"},"version":0}
%---
%[output:43f818e8]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"0.7500"},"version":0}
%---
%[output:3cedce63]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"0.8750"},"version":0}
%---
%[output:0512aa81]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"1.0000"},"version":0}
%---
%[output:8eeb2360]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"1.1250"},"version":0}
%---
%[output:9da3be31]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"1.2500"},"version":0}
%---
%[output:21d48ef2]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"1.3750"},"version":0}
%---
%[output:4909baa8]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"1.5000"},"version":0}
%---
%[output:5174914a]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"1.6250"},"version":0}
%---
%[output:6850a378]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"1.7500"},"version":0}
%---
%[output:062a608f]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"1.8750"},"version":0}
%---
%[output:93e9ac50]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"2.0000"},"version":0}
%---
%[output:2bc75bcf]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"2.1250"},"version":0}
%---
%[output:3e15590b]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"2.2500"},"version":0}
%---
%[output:5ad63e60]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"2.3750"},"version":0}
%---
%[output:2f96b81a]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"2.5000"},"version":0}
%---
%[output:4ec6be17]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"2.6250"},"version":0}
%---
%[output:4329c6a7]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"2.7500"},"version":0}
%---
%[output:1e6b9578]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"2.8750"},"version":0}
%---
%[output:8c461895]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"3.0000"},"version":0}
%---
%[output:4564161c]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"3.1250"},"version":0}
%---
%[output:399942c0]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"3.2500"},"version":0}
%---
%[output:52568e74]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"3.3750"},"version":0}
%---
%[output:9c1b4581]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"3.5000"},"version":0}
%---
%[output:23eb66a0]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"3.6250"},"version":0}
%---
%[output:88845499]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"3.7500"},"version":0}
%---
%[output:56e4ebf8]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"3.8750"},"version":0}
%---
%[output:674c21b1]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"4.0000"},"version":0}
%---
%[output:50180601]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"4.1250"},"version":0}
%---
%[output:4425a3f5]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"4.2500"},"version":0}
%---
%[output:9e0f0288]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"4.3750"},"version":0}
%---
%[output:8cfeca85]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"4.5000"},"version":0}
%---
%[output:72be3094]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"4.6250"},"version":0}
%---
%[output:9f7ff3e5]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"4.7500"},"version":0}
%---
%[output:65aff999]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"4.8750"},"version":0}
%---
%[output:66ecea3f]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"5"},"version":0}
%---
%[output:8b3e0ea9]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"5.1250"},"version":0}
%---
%[output:1e9f7840]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"5.2500"},"version":0}
%---
%[output:29c58d41]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"5.3750"},"version":0}
%---
%[output:799d7fac]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"5.5000"},"version":0}
%---
%[output:70b47185]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"5.6250"},"version":0}
%---
%[output:0a877eae]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"5.7500"},"version":0}
%---
%[output:8623ed0d]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"5.8750"},"version":0}
%---
%[output:71db15d8]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"6.0000"},"version":0}
%---
%[output:97bd5390]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"6.1250"},"version":0}
%---
%[output:4c04bdb9]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"6.2500"},"version":0}
%---
%[output:3fd60ce3]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"6.3750"},"version":0}
%---
%[output:443f10fc]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"6.5000"},"version":0}
%---
%[output:5937b00c]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"6.6250"},"version":0}
%---
%[output:890d5137]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"6.7500"},"version":0}
%---
%[output:8359d3ef]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"6.8750"},"version":0}
%---
%[output:8f02e387]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"7.0000"},"version":0}
%---
%[output:21792d8c]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"7.1250"},"version":0}
%---
%[output:139229e6]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"7.2500"},"version":0}
%---
%[output:9093d9a2]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"7.3750"},"version":0}
%---
%[output:0bf9fdad]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"7.5000"},"version":0}
%---
%[output:3c14fac7]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"7.6250"},"version":0}
%---
%[output:43fb5f33]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"7.7500"},"version":0}
%---
%[output:1b4ae426]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"7.8750"},"version":0}
%---
%[output:2c9c9ab7]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"8.0000"},"version":0}
%---
%[output:6a73c26d]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"8.1250"},"version":0}
%---
%[output:7102f20f]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"8.2500"},"version":0}
%---
%[output:88c830a0]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"8.3750"},"version":0}
%---
%[output:905ae766]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"8.5000"},"version":0}
%---
%[output:0b267d00]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"8.6250"},"version":0}
%---
%[output:1574076d]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"8.7500"},"version":0}
%---
%[output:5d22b67e]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"8.8750"},"version":0}
%---
%[output:76bd1a6c]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"9.0000"},"version":0}
%---
%[output:9d410eca]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"9.1250"},"version":0}
%---
%[output:384efeb9]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"9.2500"},"version":0}
%---
%[output:1252cc27]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"9.3750"},"version":0}
%---
%[output:76bb0a29]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"9.5000"},"version":0}
%---
%[output:1b65c340]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"9.6250"},"version":0}
%---
%[output:328d6f15]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"9.7500"},"version":0}
%---
%[output:52181767]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"9.8750"},"version":0}
%---
%[output:9f724b22]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"10"},"version":0}
%---
%[output:4cb6043d]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"10.1250"},"version":0}
%---
%[output:33c84327]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"10.2500"},"version":0}
%---
%[output:1568c1f2]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"10.3750"},"version":0}
%---
%[output:2a92733c]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"10.5000"},"version":0}
%---
%[output:45aa9910]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"10.6250"},"version":0}
%---
%[output:9e61474c]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"10.7500"},"version":0}
%---
%[output:36bb927b]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"10.8750"},"version":0}
%---
%[output:8999f1a0]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"11.0000"},"version":0}
%---
%[output:64e8eaab]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"11.1250"},"version":0}
%---
%[output:50fee51b]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"11.2500"},"version":0}
%---
%[output:063188c6]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"11.3750"},"version":0}
%---
%[output:1f9feb4f]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"11.5000"},"version":0}
%---
%[output:9ce6d470]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"11.6250"},"version":0}
%---
%[output:599f3320]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"11.7500"},"version":0}
%---
%[output:8b89ab8d]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"11.8750"},"version":0}
%---
%[output:011a88b0]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"12.0000"},"version":0}
%---
%[output:3214ef59]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"12.1250"},"version":0}
%---
%[output:6cc52d94]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"12.2500"},"version":0}
%---
%[output:68743e90]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"12.3750"},"version":0}
%---
%[output:78c8c4cd]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"12.5000"},"version":0}
%---
%[output:07d9edb2]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"12.6250"},"version":0}
%---
%[output:22cd13c0]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"12.7500"},"version":0}
%---
%[output:43fc3e51]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"12.8750"},"version":0}
%---
%[output:9c915a32]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"13.0000"},"version":0}
%---
%[output:96359e81]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"13.1250"},"version":0}
%---
%[output:38d1403a]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"13.2500"},"version":0}
%---
%[output:33658f73]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"13.3750"},"version":0}
%---
%[output:0e7778ec]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"13.5000"},"version":0}
%---
%[output:0d19b4f5]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"13.6250"},"version":0}
%---
%[output:9982c429]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"13.7500"},"version":0}
%---
%[output:689687c7]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"13.8750"},"version":0}
%---
%[output:9a4a49da]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"14.0000"},"version":0}
%---
%[output:7cfd1787]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"14.1250"},"version":0}
%---
%[output:5e63efe4]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"14.2500"},"version":0}
%---
%[output:223a66be]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"14.3750"},"version":0}
%---
%[output:6a935bea]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"14.5000"},"version":0}
%---
%[output:3371492d]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"14.6250"},"version":0}
%---
%[output:870790e4]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"14.7500"},"version":0}
%---
%[output:34538248]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"14.8750"},"version":0}
%---
%[output:3df6f7f0]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"15"},"version":0}
%---
%[output:28016d60]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"15.1250"},"version":0}
%---
%[output:1aadc492]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"15.2500"},"version":0}
%---
%[output:694a84f7]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"15.3750"},"version":0}
%---
%[output:13ce17c3]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"15.5000"},"version":0}
%---
%[output:37358051]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"15.6250"},"version":0}
%---
%[output:7af49e47]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"15.7500"},"version":0}
%---
%[output:5bda90ad]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"15.8750"},"version":0}
%---
%[output:02d3e201]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"16"},"version":0}
%---
%[output:6f1db92d]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"16.1250"},"version":0}
%---
%[output:91134c5b]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"16.2500"},"version":0}
%---
%[output:59afc586]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"16.3750"},"version":0}
%---
%[output:22f9fa0e]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"16.5000"},"version":0}
%---
%[output:9da7cd9c]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"16.6250"},"version":0}
%---
%[output:1b83bb7e]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"16.7500"},"version":0}
%---
%[output:4751edfe]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"16.8750"},"version":0}
%---
%[output:024c3c35]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"17"},"version":0}
%---
%[output:4c3eb3b0]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"17.1250"},"version":0}
%---
%[output:51756178]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"17.2500"},"version":0}
%---
%[output:4b09fc72]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"17.3750"},"version":0}
%---
%[output:8dba2d56]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"17.5000"},"version":0}
%---
%[output:139759b8]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"17.6250"},"version":0}
%---
%[output:8ccb9a5c]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"17.7500"},"version":0}
%---
%[output:99724e8f]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"17.8750"},"version":0}
%---
%[output:1ae99635]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"18.0000"},"version":0}
%---
%[output:60393c8d]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"18.1250"},"version":0}
%---
%[output:7d5728bd]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"18.2500"},"version":0}
%---
%[output:54d1675b]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"18.3750"},"version":0}
%---
%[output:069a3da4]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"18.5000"},"version":0}
%---
%[output:83ff8d41]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"18.6250"},"version":0}
%---
%[output:5d696025]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"18.7500"},"version":0}
%---
%[output:01b69dd4]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"18.8750"},"version":0}
%---
%[output:535f2d3e]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"19"},"version":0}
%---
%[output:75e40f95]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"19.1250"},"version":0}
%---
%[output:78200454]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"19.2500"},"version":0}
%---
%[output:76d7aefc]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"19.3750"},"version":0}
%---
%[output:215e545d]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"19.5000"},"version":0}
%---
%[output:2e22e3f5]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"19.6250"},"version":0}
%---
%[output:03de8895]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"19.7500"},"version":0}
%---
%[output:09cba00b]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"19.8750"},"version":0}
%---
%[output:7b65ba6e]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"20"},"version":0}
%---
%[output:77b435d0]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"20.1250"},"version":0}
%---
%[output:75272ce8]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"20.2500"},"version":0}
%---
%[output:84a92e87]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"20.3750"},"version":0}
%---
%[output:2d215bbb]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"20.5000"},"version":0}
%---
%[output:0c4442d6]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"20.6250"},"version":0}
%---
%[output:964605ac]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"20.7500"},"version":0}
%---
%[output:6aacb034]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"20.8750"},"version":0}
%---
%[output:7a8c0d75]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"21.0000"},"version":0}
%---
%[output:6555e12f]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"21.1250"},"version":0}
%---
%[output:6c1e140e]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"21.2500"},"version":0}
%---
%[output:9c9cc535]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"21.3750"},"version":0}
%---
%[output:7bd66a30]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"21.5000"},"version":0}
%---
%[output:5b9e7555]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"21.6250"},"version":0}
%---
%[output:1d98a0ea]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"21.7500"},"version":0}
%---
%[output:5fd96aae]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"21.8750"},"version":0}
%---
%[output:706ec828]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"22.0000"},"version":0}
%---
%[output:9b4de024]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"22.1250"},"version":0}
%---
%[output:2c02632f]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"22.2500"},"version":0}
%---
%[output:226614ba]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"22.3750"},"version":0}
%---
%[output:514b8679]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"22.5000"},"version":0}
%---
%[output:5a0921bd]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"22.6250"},"version":0}
%---
%[output:6eb7f962]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"22.7500"},"version":0}
%---
%[output:25f9a8f7]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"22.8750"},"version":0}
%---
%[output:1fcd2943]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"23"},"version":0}
%---
%[output:613dfbc8]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"23.1250"},"version":0}
%---
%[output:6fb2c1d3]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"23.2500"},"version":0}
%---
%[output:6a77544f]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"23.3750"},"version":0}
%---
%[output:48d54184]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"23.5000"},"version":0}
%---
%[output:1ab92c76]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"23.6250"},"version":0}
%---
%[output:2b17c194]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"23.7500"},"version":0}
%---
%[output:3d6322a2]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"23.8750"},"version":0}
%---
%[output:4eba9adf]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"24"},"version":0}
%---
%[output:01fea283]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"24.1250"},"version":0}
%---
%[output:84bda8bb]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"24.2500"},"version":0}
%---
%[output:88b358bf]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"24.3750"},"version":0}
%---
%[output:80fd3a13]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"24.5000"},"version":0}
%---
%[output:9a8f6e32]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"24.6250"},"version":0}
%---
%[output:0a4efcb4]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"24.7500"},"version":0}
%---
%[output:81bec6cd]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"24.8750"},"version":0}
%---
%[output:18be26ce]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"25"},"version":0}
%---
%[output:9656efa4]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"25.1250"},"version":0}
%---
%[output:1b67f9a1]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"25.2500"},"version":0}
%---
%[output:2bfe19c1]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"25.3750"},"version":0}
%---
%[output:8a3f7eb6]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"25.5000"},"version":0}
%---
%[output:38b7a692]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"25.6250"},"version":0}
%---
%[output:2929b9e0]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"25.7500"},"version":0}
%---
%[output:1880a2b7]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"25.8750"},"version":0}
%---
%[output:8211b003]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"26"},"version":0}
%---
%[output:1021cb98]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"26.1250"},"version":0}
%---
%[output:136342e5]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"26.2500"},"version":0}
%---
%[output:00c8a90d]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"26.3750"},"version":0}
%---
%[output:5bc653c5]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"26.5000"},"version":0}
%---
%[output:52e165b8]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"26.6250"},"version":0}
%---
%[output:69718805]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"26.7500"},"version":0}
%---
%[output:4454b93e]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"26.8750"},"version":0}
%---
%[output:2eb0fc25]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"27"},"version":0}
%---
%[output:2f1ebc73]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"27.1250"},"version":0}
%---
%[output:6ac1855c]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"27.2500"},"version":0}
%---
%[output:46048279]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"27.3750"},"version":0}
%---
%[output:5a9134e4]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"27.5000"},"version":0}
%---
%[output:75bec839]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"27.6250"},"version":0}
%---
%[output:266dd468]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"27.7500"},"version":0}
%---
%[output:7b45b5e8]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"27.8750"},"version":0}
%---
%[output:4c6b05e4]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"28.0000"},"version":0}
%---
%[output:6a47fc37]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"28.1250"},"version":0}
%---
%[output:8c1ecd7e]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"28.2500"},"version":0}
%---
%[output:363e31cc]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"28.3750"},"version":0}
%---
%[output:165731f1]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"28.5000"},"version":0}
%---
%[output:0eeb0472]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"28.6250"},"version":0}
%---
%[output:1b195a90]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"28.7500"},"version":0}
%---
%[output:270a5f25]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"28.8750"},"version":0}
%---
%[output:44414d94]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"29.0000"},"version":0}
%---
%[output:23330413]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"29.1250"},"version":0}
%---
%[output:3a6cb781]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"29.2500"},"version":0}
%---
%[output:95fbffbe]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"29.3750"},"version":0}
%---
%[output:64923a30]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"29.5000"},"version":0}
%---
%[output:99673cf8]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"29.6250"},"version":0}
%---
%[output:0c699fe3]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"29.7500"},"version":0}
%---
%[output:897f6040]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"29.8750"},"version":0}
%---
%[output:9cdbd7e5]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"30.0000"},"version":0}
%---
%[output:35624b29]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"30.1250"},"version":0}
%---
%[output:23ecd33e]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"30.2500"},"version":0}
%---
%[output:030acaa6]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"30.3750"},"version":0}
%---
%[output:54a60db0]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"30.5000"},"version":0}
%---
%[output:95e948af]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"30.6250"},"version":0}
%---
%[output:45fa1c6c]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"30.7500"},"version":0}
%---
%[output:58cff72d]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"30.8750"},"version":0}
%---
%[output:274d8d71]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"31.0000"},"version":0}
%---
%[output:66668248]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"31.1250"},"version":0}
%---
%[output:07fe4094]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"31.2500"},"version":0}
%---
%[output:5bbb967a]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"31.3750"},"version":0}
%---
%[output:4d25b75b]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"31.5000"},"version":0}
%---
%[output:603fdb6a]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"31.6250"},"version":0}
%---
%[output:955a827a]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"31.7500"},"version":0}
%---
%[output:755e79a7]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"31.8750"},"version":0}
%---
%[output:40e509b8]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"32.0000"},"version":0}
%---
%[output:321f1155]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"32.1250"},"version":0}
%---
%[output:6f92ed21]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"32.2500"},"version":0}
%---
%[output:7b672439]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"32.3750"},"version":0}
%---
%[output:1ab9803a]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"32.5000"},"version":0}
%---
%[output:149db7da]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"32.6250"},"version":0}
%---
%[output:0b5b4b5a]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"32.7500"},"version":0}
%---
%[output:1c31cd71]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"32.8750"},"version":0}
%---
%[output:2f6f7c77]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"33.0000"},"version":0}
%---
%[output:22091076]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"33.1250"},"version":0}
%---
%[output:47d5fa64]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"33.2500"},"version":0}
%---
%[output:74c79afc]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"33.3750"},"version":0}
%---
%[output:96ff3f78]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"33.5000"},"version":0}
%---
%[output:801a9c5b]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"33.6250"},"version":0}
%---
%[output:132d41c4]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"33.7500"},"version":0}
%---
%[output:81cc19d7]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"33.8750"},"version":0}
%---
%[output:1696d515]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"34.0000"},"version":0}
%---
%[output:37fb1110]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"34.1250"},"version":0}
%---
%[output:31aa9287]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"34.2500"},"version":0}
%---
%[output:937994ee]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"34.3750"},"version":0}
%---
%[output:7b9b84f7]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"34.5000"},"version":0}
%---
%[output:78338b7b]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"34.6250"},"version":0}
%---
%[output:8ea3a2c2]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"34.7500"},"version":0}
%---
%[output:8c20c088]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"34.8750"},"version":0}
%---
%[output:7fb82433]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"35.0000"},"version":0}
%---
%[output:396ebf35]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"35.1250"},"version":0}
%---
%[output:7e334680]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"35.2500"},"version":0}
%---
%[output:04c5d578]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"35.3750"},"version":0}
%---
%[output:50b9dcd4]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"35.5000"},"version":0}
%---
%[output:0c3c5267]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"35.6250"},"version":0}
%---
%[output:82408e79]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"35.7500"},"version":0}
%---
%[output:381f08aa]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"35.8750"},"version":0}
%---
%[output:65403290]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"36"},"version":0}
%---
%[output:059555df]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"36.1250"},"version":0}
%---
%[output:5b3c959b]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"36.2500"},"version":0}
%---
%[output:110c4a7c]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"36.3750"},"version":0}
%---
%[output:6fcf1ccf]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"36.5000"},"version":0}
%---
%[output:545316d8]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"36.6250"},"version":0}
%---
%[output:4e8e2520]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"36.7500"},"version":0}
%---
%[output:44fb241d]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"36.8750"},"version":0}
%---
%[output:71148cef]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"37"},"version":0}
%---
%[output:3b4a7852]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"37.1250"},"version":0}
%---
%[output:6823c738]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"37.2500"},"version":0}
%---
%[output:7a0d0ca7]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"37.3750"},"version":0}
%---
%[output:9ef39bf5]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"37.5000"},"version":0}
%---
%[output:51d8c4b6]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"37.6250"},"version":0}
%---
%[output:6192d30b]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"37.7500"},"version":0}
%---
%[output:3c7ccce6]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"37.8750"},"version":0}
%---
%[output:177d8c00]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"38"},"version":0}
%---
%[output:356bcfe5]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"38.1250"},"version":0}
%---
%[output:28080b05]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"38.2500"},"version":0}
%---
%[output:9a14a5bb]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"38.3750"},"version":0}
%---
%[output:28ade82f]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"38.5000"},"version":0}
%---
%[output:1a641902]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"38.6250"},"version":0}
%---
%[output:28593b5b]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"38.7500"},"version":0}
%---
%[output:70c22289]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"38.8750"},"version":0}
%---
%[output:967637ea]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"39.0000"},"version":0}
%---
%[output:7a96a34e]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"39.1250"},"version":0}
%---
%[output:7885d04c]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"39.2500"},"version":0}
%---
%[output:17e623ac]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"39.3750"},"version":0}
%---
%[output:9d4942fb]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"39.5000"},"version":0}
%---
%[output:46c0b5dd]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"39.6250"},"version":0}
%---
%[output:39dccf90]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"39.7500"},"version":0}
%---
%[output:8fe6fe1b]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"39.8750"},"version":0}
%---
%[output:2eda21b3]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"40.0000"},"version":0}
%---
%[output:71672212]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"40.1250"},"version":0}
%---
%[output:69273952]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"40.2500"},"version":0}
%---
%[output:3f6f2ce7]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"40.3750"},"version":0}
%---
%[output:058742a4]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"40.5000"},"version":0}
%---
%[output:431c0cb8]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"40.6250"},"version":0}
%---
%[output:53b6ab01]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"40.7500"},"version":0}
%---
%[output:27cafebb]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"40.8750"},"version":0}
%---
%[output:18da230d]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"41.0000"},"version":0}
%---
%[output:6fea10df]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"41.1250"},"version":0}
%---
%[output:96d15cf4]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"41.2500"},"version":0}
%---
%[output:53c24019]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"41.3750"},"version":0}
%---
%[output:7f702b9d]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"41.5000"},"version":0}
%---
%[output:28180a47]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"41.6250"},"version":0}
%---
%[output:35aa60ec]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"41.7500"},"version":0}
%---
%[output:9661b423]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"41.8750"},"version":0}
%---
%[output:9f1422bf]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"42"},"version":0}
%---
%[output:6754ab06]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"42.1250"},"version":0}
%---
%[output:956eeda1]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"42.2500"},"version":0}
%---
%[output:7651eee1]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"42.3750"},"version":0}
%---
%[output:92380a8b]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"42.5000"},"version":0}
%---
%[output:6be12667]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"42.6250"},"version":0}
%---
%[output:60599c40]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"42.7500"},"version":0}
%---
%[output:78e0c533]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"42.8750"},"version":0}
%---
%[output:84ddf029]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"43"},"version":0}
%---
%[output:4b178b3d]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"43.1250"},"version":0}
%---
%[output:8863680a]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"43.2500"},"version":0}
%---
%[output:1d45b848]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"43.3750"},"version":0}
%---
%[output:3639c26b]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"43.5000"},"version":0}
%---
%[output:87831e8d]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"43.6250"},"version":0}
%---
%[output:8477d437]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"43.7500"},"version":0}
%---
%[output:4a3082f5]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"43.8750"},"version":0}
%---
%[output:95e26c08]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"44"},"version":0}
%---
%[output:7004646e]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"44.1250"},"version":0}
%---
%[output:11612855]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"44.2500"},"version":0}
%---
%[output:71305314]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"44.3750"},"version":0}
%---
%[output:06f8ded7]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"44.5000"},"version":0}
%---
%[output:1a641fef]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"44.6250"},"version":0}
%---
%[output:569539f6]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"44.7500"},"version":0}
%---
%[output:264acde4]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"44.8750"},"version":0}
%---
%[output:76b6f09b]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"45"},"version":0}
%---
%[output:620bc512]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"45.1250"},"version":0}
%---
%[output:7da73d49]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"45.2500"},"version":0}
%---
%[output:1885aa20]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"45.3750"},"version":0}
%---
%[output:0c2fc8e8]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"45.5000"},"version":0}
%---
%[output:7d4eaf7c]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"45.6250"},"version":0}
%---
%[output:6563ea6b]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"45.7500"},"version":0}
%---
%[output:2a17c920]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"45.8750"},"version":0}
%---
%[output:0314a60f]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"46"},"version":0}
%---
%[output:6a43efc5]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"46.1250"},"version":0}
%---
%[output:2431247b]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"46.2500"},"version":0}
%---
%[output:1b85c110]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"46.3750"},"version":0}
%---
%[output:744a84f8]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"46.5000"},"version":0}
%---
%[output:6205f959]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"46.6250"},"version":0}
%---
%[output:8407e28e]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"46.7500"},"version":0}
%---
%[output:9ea3a873]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"46.8750"},"version":0}
%---
%[output:24765adb]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"47"},"version":0}
%---
%[output:019b01a0]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"47.1250"},"version":0}
%---
%[output:5bf46f37]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"47.2500"},"version":0}
%---
%[output:58acdc61]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"47.3750"},"version":0}
%---
%[output:78e17c4e]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"47.5000"},"version":0}
%---
%[output:26c3f164]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"47.6250"},"version":0}
%---
%[output:8257d51a]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"47.7500"},"version":0}
%---
%[output:41e6b9b2]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"47.8750"},"version":0}
%---
%[output:93f978bc]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"48.0000"},"version":0}
%---
%[output:99db324d]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"48.1250"},"version":0}
%---
%[output:838e2deb]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"48.2500"},"version":0}
%---
%[output:489ee308]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"48.3750"},"version":0}
%---
%[output:381d25ba]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"48.5000"},"version":0}
%---
%[output:3732681b]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"48.6250"},"version":0}
%---
%[output:395eedaa]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"48.7500"},"version":0}
%---
%[output:9f2efd6d]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"48.8750"},"version":0}
%---
%[output:0b064d9d]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"49"},"version":0}
%---
%[output:0c3b85e8]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"49.1250"},"version":0}
%---
%[output:89b380fd]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"49.2500"},"version":0}
%---
%[output:0fca8b88]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"49.3750"},"version":0}
%---
%[output:228e289e]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"49.5000"},"version":0}
%---
%[output:6943a535]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"49.6250"},"version":0}
%---
%[output:305a5172]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"49.7500"},"version":0}
%---
%[output:1ee3b48e]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"49.8750"},"version":0}
%---
%[output:56b344e8]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"50.0000"},"version":0}
%---
%[output:5f4d8a5b]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"50.1250"},"version":0}
%---
%[output:964b5714]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"50.2500"},"version":0}
%---
%[output:9dbbea33]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"50.3750"},"version":0}
%---
%[output:10a0b4f5]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"50.5000"},"version":0}
%---
%[output:760d7c0c]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"50.6250"},"version":0}
%---
%[output:21f6e935]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"50.7500"},"version":0}
%---
%[output:69df8682]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"50.8750"},"version":0}
%---
%[output:79e75d23]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"51"},"version":0}
%---
%[output:5b6d7efd]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"51.1250"},"version":0}
%---
%[output:85a938c2]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"51.2500"},"version":0}
%---
%[output:1ef591da]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"51.3750"},"version":0}
%---
%[output:324f5a68]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"51.5000"},"version":0}
%---
%[output:05d2a5cf]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"51.6250"},"version":0}
%---
%[output:82d6592f]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"51.7500"},"version":0}
%---
%[output:255ad6d8]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"51.8750"},"version":0}
%---
%[output:5289743b]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"52.0000"},"version":0}
%---
%[output:85401496]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"52.1250"},"version":0}
%---
%[output:40b738d3]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"52.2500"},"version":0}
%---
%[output:1bd329c7]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"52.3750"},"version":0}
%---
%[output:3763fbd1]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"52.5000"},"version":0}
%---
%[output:96ee24d6]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"52.6250"},"version":0}
%---
%[output:99e01ee0]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"52.7500"},"version":0}
%---
%[output:3bc25895]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"52.8750"},"version":0}
%---
%[output:1baae318]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"53"},"version":0}
%---
%[output:0c76f1d0]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"53.1250"},"version":0}
%---
%[output:4945bf43]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"53.2500"},"version":0}
%---
%[output:96e18e9a]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"53.3750"},"version":0}
%---
%[output:14ec4758]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"53.5000"},"version":0}
%---
%[output:69f6ec07]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"53.6250"},"version":0}
%---
%[output:9915aa92]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"53.7500"},"version":0}
%---
%[output:2060f380]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"53.8750"},"version":0}
%---
%[output:5e2902ee]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"54.0000"},"version":0}
%---
%[output:1406a359]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"54.1250"},"version":0}
%---
%[output:564aba1d]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"54.2500"},"version":0}
%---
%[output:3e3baa49]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"54.3750"},"version":0}
%---
%[output:548e9430]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"54.5000"},"version":0}
%---
%[output:7ad947db]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"54.6250"},"version":0}
%---
%[output:3b460793]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"54.7500"},"version":0}
%---
%[output:84f61f1f]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"54.8750"},"version":0}
%---
%[output:151f66f7]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"55.0000"},"version":0}
%---
%[output:9af9a3a3]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"55.1250"},"version":0}
%---
%[output:191448a2]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"55.2500"},"version":0}
%---
%[output:9db83ffc]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"55.3750"},"version":0}
%---
%[output:7127c4b6]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"55.5000"},"version":0}
%---
%[output:00306098]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"55.6250"},"version":0}
%---
%[output:0b054eb4]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"55.7500"},"version":0}
%---
%[output:0e436b1d]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"55.8750"},"version":0}
%---
%[output:01ba0877]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"56.0000"},"version":0}
%---
%[output:8ed2f7b4]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"56.1250"},"version":0}
%---
%[output:072b8c96]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"56.2500"},"version":0}
%---
%[output:417e144e]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"56.3750"},"version":0}
%---
%[output:9085762e]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"56.5000"},"version":0}
%---
%[output:483f17bc]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"56.6250"},"version":0}
%---
%[output:34efc0bc]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"56.7500"},"version":0}
%---
%[output:47c8f7d7]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"56.8750"},"version":0}
%---
%[output:9d56ccc0]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"57.0000"},"version":0}
%---
%[output:9035a2ce]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"57.1250"},"version":0}
%---
%[output:3f05b86c]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"57.2500"},"version":0}
%---
%[output:3f61d08d]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"57.3750"},"version":0}
%---
%[output:10961157]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"57.5000"},"version":0}
%---
%[output:5914a913]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"57.6250"},"version":0}
%---
%[output:5118c604]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"57.7500"},"version":0}
%---
%[output:48cbe5d4]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"57.8750"},"version":0}
%---
%[output:26829e64]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"58.0000"},"version":0}
%---
%[output:27b587b6]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"58.1250"},"version":0}
%---
%[output:380980d0]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"58.2500"},"version":0}
%---
%[output:6e262f29]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"58.3750"},"version":0}
%---
%[output:6defad62]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"58.5000"},"version":0}
%---
%[output:4cec0fec]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"58.6250"},"version":0}
%---
%[output:0b28b513]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"58.7500"},"version":0}
%---
%[output:407931e4]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"58.8750"},"version":0}
%---
%[output:0cff5c09]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"59.0000"},"version":0}
%---
%[output:2d32fbf1]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"59.1250"},"version":0}
%---
%[output:33b53d72]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"59.2500"},"version":0}
%---
%[output:6e4430e8]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"59.3750"},"version":0}
%---
%[output:2fac5439]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"59.5000"},"version":0}
%---
%[output:463cf34f]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"59.6250"},"version":0}
%---
%[output:77639cee]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"59.7500"},"version":0}
%---
%[output:736f75a4]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"59.8750"},"version":0}
%---
%[output:37c8c8ff]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"60"},"version":0}
%---
%[output:8fcdcf3d]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"60.1250"},"version":0}
%---
%[output:0cbb6239]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"60.2500"},"version":0}
%---
%[output:70da4a1f]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"60.3750"},"version":0}
%---
%[output:8c72b683]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"60.5000"},"version":0}
%---
%[output:7732bf28]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"60.6250"},"version":0}
%---
%[output:38309d76]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"60.7500"},"version":0}
%---
%[output:121ab693]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"60.8750"},"version":0}
%---
%[output:72479b01]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"61.0000"},"version":0}
%---
%[output:144114a2]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"61.1250"},"version":0}
%---
%[output:89768699]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"61.2500"},"version":0}
%---
%[output:859e995f]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"61.3750"},"version":0}
%---
%[output:5340c61d]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"61.5000"},"version":0}
%---
%[output:259b9533]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"61.6250"},"version":0}
%---
%[output:4d582263]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"61.7500"},"version":0}
%---
%[output:12a6e5ac]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"61.8750"},"version":0}
%---
%[output:5b33f8ce]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"62"},"version":0}
%---
%[output:804edcbe]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"62.1250"},"version":0}
%---
%[output:04cd8b18]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"62.2500"},"version":0}
%---
%[output:8ac5e673]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"62.3750"},"version":0}
%---
%[output:7f2f9789]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"62.5000"},"version":0}
%---
%[output:78836302]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"62.6250"},"version":0}
%---
%[output:1d244cce]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"62.7500"},"version":0}
%---
%[output:3cae733e]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"62.8750"},"version":0}
%---
%[output:5a98d013]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"63.0000"},"version":0}
%---
%[output:4a2f020b]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"63.1250"},"version":0}
%---
%[output:3ac4683b]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"63.2500"},"version":0}
%---
%[output:5e401a06]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"63.3750"},"version":0}
%---
%[output:1cc3633b]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"63.5000"},"version":0}
%---
%[output:970ad9af]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"63.6250"},"version":0}
%---
%[output:1950352c]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"63.7500"},"version":0}
%---
%[output:3cddb129]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"63.8750"},"version":0}
%---
%[output:120fe386]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"64"},"version":0}
%---
%[output:05c29855]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"64.1250"},"version":0}
%---
%[output:56c1b452]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"64.2500"},"version":0}
%---
%[output:20072cf5]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"64.3750"},"version":0}
%---
%[output:45194e1b]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"64.5000"},"version":0}
%---
%[output:19b62263]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"64.6250"},"version":0}
%---
%[output:994f6dfe]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"64.7500"},"version":0}
%---
%[output:426eea91]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"64.8750"},"version":0}
%---
%[output:47eb7e26]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"65.0000"},"version":0}
%---
%[output:4dba9289]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"65.1250"},"version":0}
%---
%[output:7b54bd6e]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"65.2500"},"version":0}
%---
%[output:32fa3add]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"65.3750"},"version":0}
%---
%[output:7046e0a8]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"65.5000"},"version":0}
%---
%[output:8dea4457]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"65.6250"},"version":0}
%---
%[output:18c8dcef]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"65.7500"},"version":0}
%---
%[output:9da87ee6]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"65.8750"},"version":0}
%---
%[output:6562bce7]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"66.0000"},"version":0}
%---
%[output:729b5839]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"66.1250"},"version":0}
%---
%[output:51d0d79b]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"66.2500"},"version":0}
%---
%[output:84a8bada]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"66.3750"},"version":0}
%---
%[output:94ced4b9]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"66.5000"},"version":0}
%---
%[output:4a908277]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"66.6250"},"version":0}
%---
%[output:533db7d4]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"66.7500"},"version":0}
%---
%[output:7197219b]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"66.8750"},"version":0}
%---
%[output:0314c648]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"67"},"version":0}
%---
%[output:7e27bbaa]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"67.1250"},"version":0}
%---
%[output:3843cfd1]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"67.2500"},"version":0}
%---
%[output:01e79ea5]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"67.3750"},"version":0}
%---
%[output:0ae75e26]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"67.5000"},"version":0}
%---
%[output:5ceb97da]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"67.6250"},"version":0}
%---
%[output:2852c5ef]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"67.7500"},"version":0}
%---
%[output:9511a384]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"67.8750"},"version":0}
%---
%[output:5cbb6b72]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"68.0000"},"version":0}
%---
%[output:20bb45ee]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"68.1250"},"version":0}
%---
%[output:79b6ff1a]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"68.2500"},"version":0}
%---
%[output:36dc67ec]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"68.3750"},"version":0}
%---
%[output:691d1057]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"68.5000"},"version":0}
%---
%[output:36de0900]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"68.6250"},"version":0}
%---
%[output:6098be5d]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"68.7500"},"version":0}
%---
%[output:6759a4b1]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"68.8750"},"version":0}
%---
%[output:4056dafc]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"69"},"version":0}
%---
%[output:056303c7]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"69.1250"},"version":0}
%---
%[output:49cd1bbb]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"69.2500"},"version":0}
%---
%[output:9995b870]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"69.3750"},"version":0}
%---
%[output:0e356ab9]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"69.5000"},"version":0}
%---
%[output:8f00889b]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"69.6250"},"version":0}
%---
%[output:74e19455]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"69.7500"},"version":0}
%---
%[output:8ec79a8c]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"69.8750"},"version":0}
%---
%[output:28fea522]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"70.0000"},"version":0}
%---
%[output:365f4b1e]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"70.1250"},"version":0}
%---
%[output:3d47d553]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"70.2500"},"version":0}
%---
%[output:98d7f80d]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"70.3750"},"version":0}
%---
%[output:5ab8911d]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"70.5000"},"version":0}
%---
%[output:26d022a8]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"70.6250"},"version":0}
%---
%[output:6e63325a]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"70.7500"},"version":0}
%---
%[output:9bd48408]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"70.8750"},"version":0}
%---
%[output:7d06c9f5]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"71"},"version":0}
%---
%[output:63df2a33]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"71.1250"},"version":0}
%---
%[output:6d414100]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"71.2500"},"version":0}
%---
%[output:1c072ab4]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"71.3750"},"version":0}
%---
%[output:229ec812]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"71.5000"},"version":0}
%---
%[output:41a392e7]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"71.6250"},"version":0}
%---
%[output:178c0d9d]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"71.7500"},"version":0}
%---
%[output:06d54718]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"71.8750"},"version":0}
%---
%[output:2f6077ad]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"72.0000"},"version":0}
%---
%[output:182a79a5]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"72.1250"},"version":0}
%---
%[output:22076f86]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"72.2500"},"version":0}
%---
%[output:5b538180]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"72.3750"},"version":0}
%---
%[output:78a90a2a]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"72.5000"},"version":0}
%---
%[output:8e92e969]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"72.6250"},"version":0}
%---
%[output:1d3e6eb6]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"72.7500"},"version":0}
%---
%[output:3e215f15]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"72.8750"},"version":0}
%---
%[output:2fbbf178]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"73"},"version":0}
%---
%[output:0fc2bba4]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"73.1250"},"version":0}
%---
%[output:025a2069]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"73.2500"},"version":0}
%---
%[output:95010784]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"73.3750"},"version":0}
%---
%[output:7467d443]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"73.5000"},"version":0}
%---
%[output:25231ef8]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"73.6250"},"version":0}
%---
%[output:1101d565]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"73.7500"},"version":0}
%---
%[output:778212f0]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"73.8750"},"version":0}
%---
%[output:5e610b10]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"74.0000"},"version":0}
%---
%[output:0722b6de]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"74.1250"},"version":0}
%---
%[output:4861de90]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"74.2500"},"version":0}
%---
%[output:14f02c7c]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"74.3750"},"version":0}
%---
%[output:7d88f410]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"74.5000"},"version":0}
%---
%[output:586ac50b]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"74.6250"},"version":0}
%---
%[output:6f64eba1]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"74.7500"},"version":0}
%---
%[output:80223990]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"74.8750"},"version":0}
%---
%[output:13712d29]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"75"},"version":0}
%---
%[output:0ec1931e]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"75.1250"},"version":0}
%---
%[output:81dd539c]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"75.2500"},"version":0}
%---
%[output:9fac415c]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"75.3750"},"version":0}
%---
%[output:5f116e0e]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"75.5000"},"version":0}
%---
%[output:1dc8fc8b]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"75.6250"},"version":0}
%---
%[output:23e2e1b0]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"75.7500"},"version":0}
%---
%[output:1f609d46]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"75.8750"},"version":0}
%---
%[output:673028df]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"76"},"version":0}
%---
%[output:9d2c4fa0]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"76.1250"},"version":0}
%---
%[output:178db462]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"76.2500"},"version":0}
%---
%[output:9a6c55d8]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"76.3750"},"version":0}
%---
%[output:2e105267]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"76.5000"},"version":0}
%---
%[output:8cff8852]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"76.6250"},"version":0}
%---
%[output:421cc68b]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"76.7500"},"version":0}
%---
%[output:2a0d85a2]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"76.8750"},"version":0}
%---
%[output:3353f2fd]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"77.0000"},"version":0}
%---
%[output:664c92a7]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"77.1250"},"version":0}
%---
%[output:7e658694]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"77.2500"},"version":0}
%---
%[output:15162885]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"77.3750"},"version":0}
%---
%[output:51ad56c2]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"77.5000"},"version":0}
%---
%[output:55c1b230]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"77.6250"},"version":0}
%---
%[output:81da06be]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"77.7500"},"version":0}
%---
%[output:957fabd9]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"77.8750"},"version":0}
%---
%[output:3cd8fc4b]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"78"},"version":0}
%---
%[output:9e23b446]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"78.1250"},"version":0}
%---
%[output:361f5440]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"78.2500"},"version":0}
%---
%[output:2215f5c9]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"78.3750"},"version":0}
%---
%[output:2ff88c13]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"78.5000"},"version":0}
%---
%[output:293ebc1d]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"78.6250"},"version":0}
%---
%[output:5438b0db]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"78.7500"},"version":0}
%---
%[output:8b2ee12b]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"78.8750"},"version":0}
%---
%[output:9714b09f]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"79.0000"},"version":0}
%---
%[output:424feca3]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"79.1250"},"version":0}
%---
%[output:21834e3e]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"79.2500"},"version":0}
%---
%[output:5b03f33d]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"79.3750"},"version":0}
%---
%[output:341272de]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"79.5000"},"version":0}
%---
%[output:01292d95]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"79.6250"},"version":0}
%---
%[output:8486e9a9]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"79.7500"},"version":0}
%---
%[output:51c99304]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"79.8750"},"version":0}
%---
%[output:9e2d8635]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"80"},"version":0}
%---
%[output:2657df05]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"80.1250"},"version":0}
%---
%[output:2f8562f6]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"80.2500"},"version":0}
%---
%[output:239b917a]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"80.3750"},"version":0}
%---
%[output:3830e1a1]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"80.5000"},"version":0}
%---
%[output:2d8d1fa7]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"80.6250"},"version":0}
%---
%[output:0693e091]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"80.7500"},"version":0}
%---
%[output:964cc592]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"80.8750"},"version":0}
%---
%[output:752a08cc]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"81.0000"},"version":0}
%---
%[output:8c6c6a54]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"81.1250"},"version":0}
%---
%[output:614c5ccf]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"81.2500"},"version":0}
%---
%[output:08907412]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"81.3750"},"version":0}
%---
%[output:6a20dbe7]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"81.5000"},"version":0}
%---
%[output:68113b07]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"81.6250"},"version":0}
%---
%[output:2eef1172]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"81.7500"},"version":0}
%---
%[output:56d33e48]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"81.8750"},"version":0}
%---
%[output:69be270b]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"82"},"version":0}
%---
%[output:5d3a4ddc]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"82.1250"},"version":0}
%---
%[output:83a85315]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"82.2500"},"version":0}
%---
%[output:7e321cb0]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"82.3750"},"version":0}
%---
%[output:0995325f]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"82.5000"},"version":0}
%---
%[output:9df4b145]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"82.6250"},"version":0}
%---
%[output:2b494034]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"82.7500"},"version":0}
%---
%[output:4c5520cd]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"82.8750"},"version":0}
%---
%[output:96ca4210]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"83.0000"},"version":0}
%---
%[output:2c03b7e6]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"83.1250"},"version":0}
%---
%[output:8d0b49d5]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"83.2500"},"version":0}
%---
%[output:7d8b3732]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"83.3750"},"version":0}
%---
%[output:43159cbb]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"83.5000"},"version":0}
%---
%[output:7e55382e]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"83.6250"},"version":0}
%---
%[output:12eb3bd4]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"83.7500"},"version":0}
%---
%[output:66b71042]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"83.8750"},"version":0}
%---
%[output:05c7e8e9]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"84.0000"},"version":0}
%---
%[output:98c9f9f4]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"84.1250"},"version":0}
%---
%[output:8fc83dc5]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"84.2500"},"version":0}
%---
%[output:62923a04]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"84.3750"},"version":0}
%---
%[output:07e40f0f]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"84.5000"},"version":0}
%---
%[output:73cc72ce]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"84.6250"},"version":0}
%---
%[output:16ae563c]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"84.7500"},"version":0}
%---
%[output:7ece0681]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"84.8750"},"version":0}
%---
%[output:21960f4f]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"85.0000"},"version":0}
%---
%[output:8e0387f5]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"85.1250"},"version":0}
%---
%[output:464479c2]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"85.2500"},"version":0}
%---
%[output:047489c9]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"85.3750"},"version":0}
%---
%[output:94f83b2c]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"85.5000"},"version":0}
%---
%[output:181206d1]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"85.6250"},"version":0}
%---
%[output:0100eb18]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"85.7500"},"version":0}
%---
%[output:05ab07ae]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"85.8750"},"version":0}
%---
%[output:036d51ec]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"86.0000"},"version":0}
%---
%[output:5a4d8a80]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"86.1250"},"version":0}
%---
%[output:7511b3b3]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"86.2500"},"version":0}
%---
%[output:3820b7f8]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"86.3750"},"version":0}
%---
%[output:29753cd6]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"86.5000"},"version":0}
%---
%[output:80166951]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"86.6250"},"version":0}
%---
%[output:6de4c5de]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"86.7500"},"version":0}
%---
%[output:2db7f1f8]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"86.8750"},"version":0}
%---
%[output:68c14a0b]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"87.0000"},"version":0}
%---
%[output:3fdf370d]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"87.1250"},"version":0}
%---
%[output:7e55fb39]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"87.2500"},"version":0}
%---
%[output:8c865a76]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"87.3750"},"version":0}
%---
%[output:15fe009f]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"87.5000"},"version":0}
%---
%[output:874066da]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"87.6250"},"version":0}
%---
%[output:94230be8]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"87.7500"},"version":0}
%---
%[output:8152fc03]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"87.8750"},"version":0}
%---
%[output:4ab10d19]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"88.0000"},"version":0}
%---
%[output:4d24f2f4]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"88.1250"},"version":0}
%---
%[output:20699bd6]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"88.2500"},"version":0}
%---
%[output:1fac4615]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"88.3750"},"version":0}
%---
%[output:078bc257]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"88.5000"},"version":0}
%---
%[output:087401ef]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"88.6250"},"version":0}
%---
%[output:90aebbfb]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"88.7500"},"version":0}
%---
%[output:8e9deb37]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"88.8750"},"version":0}
%---
%[output:3d474907]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"89"},"version":0}
%---
%[output:2f4a9930]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"89.1250"},"version":0}
%---
%[output:6ca79501]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"89.2500"},"version":0}
%---
%[output:7702b830]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"89.3750"},"version":0}
%---
%[output:97a66719]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"89.5000"},"version":0}
%---
%[output:0d0b52fd]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"89.6250"},"version":0}
%---
%[output:228b9af6]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"89.7500"},"version":0}
%---
%[output:14bb6bed]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"89.8750"},"version":0}
%---
%[output:1e94472b]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"90.0000"},"version":0}
%---
%[output:257b54dd]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"90.1250"},"version":0}
%---
%[output:0f3283fb]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"90.2500"},"version":0}
%---
%[output:870188e8]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"90.3750"},"version":0}
%---
%[output:393d000d]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"90.5000"},"version":0}
%---
%[output:6a0f3061]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"90.6250"},"version":0}
%---
%[output:702b1eef]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"90.7500"},"version":0}
%---
%[output:0fb12f55]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"90.8750"},"version":0}
%---
%[output:7748a852]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"91"},"version":0}
%---
%[output:19f568c0]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"91.1250"},"version":0}
%---
%[output:3a60f553]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"91.2500"},"version":0}
%---
%[output:553f48f4]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"91.3750"},"version":0}
%---
%[output:2b4d7004]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"91.5000"},"version":0}
%---
%[output:45a2957e]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"91.6250"},"version":0}
%---
%[output:621a2a92]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"91.7500"},"version":0}
%---
%[output:42454fd2]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"91.8750"},"version":0}
%---
%[output:68e11212]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"92"},"version":0}
%---
%[output:5692ca9b]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"92.1250"},"version":0}
%---
%[output:85e293fc]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"92.2500"},"version":0}
%---
%[output:72da5d47]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"92.3750"},"version":0}
%---
%[output:746fbf2e]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"92.5000"},"version":0}
%---
%[output:425da69e]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"92.6250"},"version":0}
%---
%[output:44467add]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"92.7500"},"version":0}
%---
%[output:2eb9db58]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"92.8750"},"version":0}
%---
%[output:638663c1]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"93"},"version":0}
%---
%[output:959dce23]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"93.1250"},"version":0}
%---
%[output:98df5fbd]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"93.2500"},"version":0}
%---
%[output:2facbc44]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"93.3750"},"version":0}
%---
%[output:447d832e]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"93.5000"},"version":0}
%---
%[output:5f4bb0cb]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"93.6250"},"version":0}
%---
%[output:0dbdb5cb]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"93.7500"},"version":0}
%---
%[output:3d0925f0]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"93.8750"},"version":0}
%---
%[output:48702547]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"94"},"version":0}
%---
%[output:3d0f8974]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"94.1250"},"version":0}
%---
%[output:8ca5d660]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"94.2500"},"version":0}
%---
%[output:1e9e16c4]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"94.3750"},"version":0}
%---
%[output:7f057364]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"94.5000"},"version":0}
%---
%[output:5f66681d]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"94.6250"},"version":0}
%---
%[output:0b5835b3]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"94.7500"},"version":0}
%---
%[output:30f3acd8]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"94.8750"},"version":0}
%---
%[output:6e29d86c]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"95"},"version":0}
%---
%[output:6ef18a41]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"95.1250"},"version":0}
%---
%[output:322bf010]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"95.2500"},"version":0}
%---
%[output:6c78e3aa]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"95.3750"},"version":0}
%---
%[output:68defd27]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"95.5000"},"version":0}
%---
%[output:98b53265]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"95.6250"},"version":0}
%---
%[output:09921628]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"95.7500"},"version":0}
%---
%[output:9d0ca7d6]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"95.8750"},"version":0}
%---
%[output:607e7357]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"96"},"version":0}
%---
%[output:0f9687b7]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"96.1250"},"version":0}
%---
%[output:3afab700]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"96.2500"},"version":0}
%---
%[output:1d0a6ebf]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"96.3750"},"version":0}
%---
%[output:8a259384]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"96.5000"},"version":0}
%---
%[output:507b6d0e]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"96.6250"},"version":0}
%---
%[output:35f56edf]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"96.7500"},"version":0}
%---
%[output:6f7ece41]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"96.8750"},"version":0}
%---
%[output:0ba0aa8e]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"97.0000"},"version":0}
%---
%[output:382e5f72]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"97.1250"},"version":0}
%---
%[output:669b5316]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"97.2500"},"version":0}
%---
%[output:52557c3b]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"97.3750"},"version":0}
%---
%[output:0b941c1e]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"97.5000"},"version":0}
%---
%[output:1def1396]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"97.6250"},"version":0}
%---
%[output:0d5ec1f4]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"97.7500"},"version":0}
%---
%[output:37f2f837]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"97.8750"},"version":0}
%---
%[output:50a9f906]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"98"},"version":0}
%---
%[output:1eddf2b9]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"98.1250"},"version":0}
%---
%[output:4647a7a7]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"98.2500"},"version":0}
%---
%[output:4ff59d3b]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"98.3750"},"version":0}
%---
%[output:655da3ac]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"98.5000"},"version":0}
%---
%[output:478093d0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Some output might be missing due to a network interruption. To get the missing output, rerun the script."}}
%---
%[output:4151fff2]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"99.2500"},"version":0}
%---
%[output:365b3483]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"99.3750"},"version":0}
%---
%[output:43a15dcf]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"99.5000"},"version":0}
%---
%[output:0b62e4d8]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"99.6250"},"version":0}
%---
%[output:22b01d39]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"99.7500"},"version":0}
%---
%[output:2ecf2c42]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"99.8750"},"version":0}
%---
%[output:2925a644]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"ans","rows":"1","value":"100"},"version":0}
%---
