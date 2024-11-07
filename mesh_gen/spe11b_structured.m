%% Matlab mesh
%% spe11b_structured, Created by Gmsh
%% ASCII
clear msh;
msh.nbNod = 242;
msh.POS = [
0 360 0;
840 360 0;
7560 360 0;
8400 360 0;
0 480 0;
840 480 0;
1680 480 0;
2520 480 0;
3360 480 0;
4200 480 0;
5040 480 0;
5880 480 0;
6720 480 0;
7560 480 0;
8400 480 0;
1680 600 0;
2520 600 0;
3360 600 0;
4200 600 0;
5040 600 0;
5880 600 0;
6720 600 0;
7560 600 0;
0 1080 0;
840 1080 0;
1680 1080 0;
2520 1080 0;
3360 1080 0;
4200 1080 0;
5040 1080 0;
5880 1080 0;
6720 1080 0;
7560 1080 0;
8400 1080 0;
0 1200 0;
840 1200 0;
1680 1200 0;
2520 1200 0;
3360 1200 0;
4200 1200 0;
5040 1200 0;
5880 1200 0;
6720 1200 0;
7560 1200 0;
8400 1200 0;
0 360 1;
840 360 1;
7560 360 1;
8400 360 1;
0 480 1;
840 480 1;
1680 480 1;
2520 480 1;
3360 480 1;
4200 480 1;
5040 480 1;
5880 480 1;
6720 480 1;
7560 480 1;
8400 480 1;
1680 600 1;
2520 600 1;
3360 600 1;
4200 600 1;
5040 600 1;
5880 600 1;
6720 600 1;
7560 600 1;
0 1080 1;
840 1080 1;
1680 1080 1;
2520 1080 1;
3360 1080 1;
4200 1080 1;
5040 1080 1;
5880 1080 1;
6720 1080 1;
7560 1080 1;
8400 1080 1;
0 1200 1;
840 1200 1;
1680 1200 1;
2520 1200 1;
3360 1200 1;
4200 1200 1;
5040 1200 1;
5880 1200 1;
6720 1200 1;
7560 1200 1;
8400 1200 1;
0 960 0;
840 960 0;
1680 960 0;
2520 960 0;
3360 960 0;
4200 960 0;
5040 960 0;
5880 960 0;
6720 960 0;
7560 960 0;
8400 960 0;
0 960 1;
840 960 1;
1680 960 1;
2520 960 1;
3360 960 1;
4200 960 1;
5040 960 1;
5880 960 1;
6720 960 1;
7560 960 1;
8400 960 1;
1680 360 0;
0 600 0;
840 600 0;
0 720 0;
840 720 0;
1680 720 0;
7560 840 0;
8400 840 0;
1680 360 1;
0 600 1;
840 600 1;
0 720 1;
840 720 1;
1680 720 1;
7560 840 1;
8400 840 1;
8400 600 0;
2520 720 0;
3360 720 0;
4200 720 0;
5040 720 0;
5880 720 0;
6720 720 0;
7560 720 0;
8400 720 0;
0 840 0;
840 840 0;
1680 840 0;
2520 840 0;
3360 840 0;
4200 840 0;
5040 840 0;
5880 840 0;
6720 840 0;
8400 600 1;
2520 720 1;
3360 720 1;
4200 720 1;
5040 720 1;
5880 720 1;
6720 720 1;
7560 720 1;
8400 720 1;
0 840 1;
840 840 1;
1680 840 1;
2520 840 1;
3360 840 1;
4200 840 1;
5040 840 1;
5880 840 1;
6720 840 1;
2520 0 0;
3360 0 0;
7560 0 0;
8400 0 0;
0 120 0;
840 120 0;
1680 120 0;
2520 120 0;
3360 120 0;
4200 120 0;
5040 120 0;
5880 120 0;
6720 120 0;
7560 120 0;
8400 120 0;
0 240 0;
840 240 0;
1680 240 0;
2520 240 0;
3360 240 0;
4200 240 0;
5040 240 0;
5880 240 0;
6720 240 0;
7560 240 0;
8400 240 0;
2520 360 0;
3360 360 0;
4200 360 0;
5040 360 0;
5880 360 0;
6720 360 0;
2520 0 1;
3360 0 1;
7560 0 1;
8400 0 1;
0 120 1;
840 120 1;
1680 120 1;
2520 120 1;
3360 120 1;
4200 120 1;
5040 120 1;
5880 120 1;
6720 120 1;
7560 120 1;
8400 120 1;
0 240 1;
840 240 1;
1680 240 1;
2520 240 1;
3360 240 1;
4200 240 1;
5040 240 1;
5880 240 1;
6720 240 1;
7560 240 1;
8400 240 1;
2520 360 1;
3360 360 1;
4200 360 1;
5040 360 1;
5880 360 1;
6720 360 1;
0 0 0;
840 0 0;
1680 0 0;
4200 0 0;
5040 0 0;
5880 0 0;
6720 0 0;
0 0 1;
840 0 1;
1680 0 1;
4200 0 1;
5040 0 1;
5880 0 1;
6720 0 1;
];
msh.MAX = max(msh.POS);
msh.MIN = min(msh.POS);
msh.HEXAS =[
 1 2 6 5 46 47 51 50 1
 3 4 15 14 48 49 60 59 1
 7 8 17 16 52 53 62 61 1
 8 9 18 17 53 54 63 62 1
 9 10 19 18 54 55 64 63 1
 10 11 20 19 55 56 65 64 1
 11 12 21 20 56 57 66 65 1
 12 13 22 21 57 58 67 66 1
 13 14 23 22 58 59 68 67 1
 24 25 36 35 69 70 81 80 1
 25 26 37 36 70 71 82 81 1
 26 27 38 37 71 72 83 82 1
 27 28 39 38 72 73 84 83 1
 28 29 40 39 73 74 85 84 1
 29 30 41 40 74 75 86 85 1
 30 31 42 41 75 76 87 86 1
 31 32 43 42 76 77 88 87 1
 32 33 44 43 77 78 89 88 1
 33 34 45 44 78 79 90 89 1
 91 92 25 24 102 103 70 69 2
 92 93 26 25 103 104 71 70 2
 93 94 27 26 104 105 72 71 2
 94 95 28 27 105 106 73 72 2
 95 96 29 28 106 107 74 73 2
 96 97 30 29 107 108 75 74 2
 97 98 31 30 108 109 76 75 2
 98 99 32 31 109 110 77 76 2
 99 100 33 32 110 111 78 77 2
 100 101 34 33 111 112 79 78 2
 2 113 7 6 47 121 52 51 3
 114 115 117 116 122 123 125 124 3
 115 16 118 117 123 61 126 125 3
 119 120 101 100 127 128 112 111 3
 5 6 115 114 50 51 123 122 4
 14 15 129 23 59 60 147 68 4
 16 17 130 118 61 62 148 126 4
 17 18 131 130 62 63 149 148 4
 18 19 132 131 63 64 150 149 4
 19 20 133 132 64 65 151 150 4
 20 21 134 133 65 66 152 151 4
 21 22 135 134 66 67 153 152 4
 22 23 136 135 67 68 154 153 4
 131 132 143 142 149 150 161 160 4
 136 137 120 119 154 155 128 127 4
 138 139 92 91 156 157 103 102 4
 139 140 93 92 157 158 104 103 4
 140 141 94 93 158 159 105 104 4
 141 142 95 94 159 160 106 105 4
 143 144 97 96 161 162 108 107 4
 144 145 98 97 162 163 109 108 4
 145 146 99 98 163 164 110 109 4
 146 119 100 99 164 127 111 110 4
 165 166 173 172 197 198 205 204 5
 167 168 179 178 199 200 211 210 5
 169 170 181 180 201 202 213 212 5
 170 171 182 181 202 203 214 213 5
 171 172 183 182 203 204 215 214 5
 172 173 184 183 204 205 216 215 5
 173 174 185 184 205 206 217 216 5
 174 175 186 185 206 207 218 217 5
 175 176 187 186 207 208 219 218 5
 176 177 188 187 208 209 220 219 5
 177 178 189 188 209 210 221 220 5
 178 179 190 189 210 211 222 221 5
 180 181 2 1 212 213 47 46 5
 181 182 113 2 213 214 121 47 5
 182 183 191 113 214 215 223 121 5
 183 184 192 191 215 216 224 223 5
 184 185 193 192 216 217 225 224 5
 185 186 194 193 217 218 226 225 5
 186 187 195 194 218 219 227 226 5
 187 188 196 195 219 220 228 227 5
 188 189 3 196 220 221 48 228 5
 189 190 4 3 221 222 49 48 5
 113 191 8 7 121 223 53 52 5
 191 192 9 8 223 224 54 53 5
 192 193 10 9 224 225 55 54 5
 193 194 11 10 225 226 56 55 5
 194 195 12 11 226 227 57 56 5
 195 196 13 12 227 228 58 57 5
 196 3 14 13 228 48 59 58 5
 23 129 137 136 68 147 155 154 5
 116 117 139 138 124 125 157 156 5
 117 118 140 139 125 126 158 157 5
 118 130 141 140 126 148 159 158 5
 130 131 142 141 148 149 160 159 5
 132 133 144 143 150 151 162 161 5
 133 134 145 144 151 152 163 162 5
 134 135 146 145 152 153 164 163 5
 135 136 119 146 153 154 127 164 5
 6 7 16 115 51 52 61 123 6
 142 143 96 95 160 161 107 106 6
 229 230 170 169 236 237 202 201 7
 230 231 171 170 237 238 203 202 7
 231 165 172 171 238 197 204 203 7
 166 232 174 173 198 239 206 205 7
 232 233 175 174 239 240 207 206 7
 233 234 176 175 240 241 208 207 7
 234 235 177 176 241 242 209 208 7
 235 167 178 177 242 199 210 209 7
];
