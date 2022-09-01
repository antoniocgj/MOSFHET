#include "mosfhet.h"

#ifdef TORUS32
#define TORUS_SCALE 4294967295.0
#else
#define TORUS_SCALE 18446744073709551615.0
#endif

double torus2double(Torus x){
  return ((double) x)/TORUS_SCALE;
}

Torus double2torus(double x){
  return (Torus) ((int64_t) (((double)TORUS_SCALE)*x));
}

/* return round(in[i] * 2^log_scale) */
uint64_t torus2int(Torus x, int log_scale){
  const uint64_t bit_size = sizeof(Torus) * 8;
  const Torus round_offset = 1UL << (bit_size - log_scale - 1);
  return (x + round_offset)>>(bit_size - log_scale);
}

/* return round(in[i] / 2^log_scale) */
Torus int2torus(uint64_t x, int log_scale){
  const uint64_t bit_size = sizeof(Torus) * 8;
  return x << (bit_size - log_scale);
}

// Random generation

#ifndef PORTABLE_BUILD
// TODO: add code src.
void generate_rnd_seed(uint64_t * p){
  if(0 == _rdrand64_step ((unsigned long long *) p) ||
     0 == _rdrand64_step ((unsigned long long *) &(p[1])) ||
     0 == _rdrand64_step ((unsigned long long *) &(p[2])) ||
     0 == _rdrand64_step ((unsigned long long *) &(p[3]))){
    printf("Random Generation Failed\n");
    return;
  }
}
#else 
void generate_rnd_seed(uint64_t * p){
  FILE *fp;
  fp = fopen("/dev/urandom", "r");
  fread(p, 1, 32, fp);
  fclose(fp);
}
#endif


#include "sha3/fips202.h"

void get_rnd_from_hash(uint64_t amount, uint8_t * pointer){
  uint64_t rnd[4];
  generate_rnd_seed(rnd);
  shake256(pointer, amount, (uint8_t *) rnd, 32);
}

void get_rnd_from_buffer(uint64_t amount, uint8_t * pointer){
  static uint8_t buffer[1024];
  static int idx = 1024;
  if(amount > (1024 - idx)){
    idx = 0;
    get_rnd_from_hash(1024, buffer);
  }
  memcpy(pointer, buffer+idx, amount);
  idx += amount;
}

void generate_random_bytes(uint64_t amount, uint8_t * pointer){
  if(amount < 512) get_rnd_from_buffer(amount, pointer);
  else get_rnd_from_hash(amount, pointer);
}

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif
double generate_normal_random(double sigma){
  uint64_t rnd[2];
  generate_random_bytes(16, (uint8_t *) rnd);
  return cos(2.*M_PI*torus2double(rnd[0]))*sqrt(-2.*log(torus2double(rnd[1])))*sigma;
}

void generate_torus_normal_random_array(Torus * out, double sigma, int N){
  for (size_t i = 0; i < N; i++){
    out[i] = double2torus(generate_normal_random(sigma));
  }
}

// Mem alloc

uint64_t _glb_mem_count = 0;
// safe_malloc
// https://stackoverflow.com/questions/48043811/creating-a-function-to-check-if-malloc-succeeded
void * safe_malloc(size_t size){
  void *ptr = malloc(size);
  if (!ptr && (size > 0)) {
    perror("malloc failed!");
    exit(EXIT_FAILURE);
  }
  // memset(ptr, 0, size); 
  return ptr;
}

void * safe_aligned_malloc(size_t size){
  void * ptr;
  #ifdef AVX512_OPT
  int err = posix_memalign(&ptr, 64, size);
  #else
  int err = posix_memalign(&ptr, 32, size);
  #endif
  if (err || (!ptr && (size > 0))) {
    perror("aligned malloc failed!");
    exit(EXIT_FAILURE);
  }
  // memset(ptr, 0, size);
  return ptr;
}

uint64_t inverse_mod(uint64_t input, uint64_t modulus) ;

uint16_t * pre_compute_odd_inverse_mod(uint16_t modulus){
  uint16_t * res = (uint16_t *) safe_aligned_malloc(sizeof(uint16_t)*modulus/2);
  for (size_t i = 0; i < modulus; i++){
    if(!(i&1)) continue;
    res[i/2] = inverse_mod(i, modulus); 
  }
  return res;
}

// Return the inverse of x mod N, for odd x and N in {512, 1024, 2048}
uint16_t inverse_mod_2N(uint16_t x, uint16_t N){
  static uint16_t N256[256] = {1, 171, 205, 439, 57, 419, 197, 239, 241, 27, 317, 423, 41, 19, 53, 479, 481, 395, 429, 407, 25, 131, 421, 207, 209, 251, 29, 391, 9, 243, 277, 447, 449, 107, 141, 375, 505, 355, 133, 175, 177, 475, 253, 359, 489, 467, 501, 415, 417, 331, 365, 343, 473, 67, 357, 143, 145, 187, 477, 327, 457, 179, 213, 383, 385, 43, 77, 311, 441, 291, 69, 111, 113, 411, 189, 295, 425, 403, 437, 351, 353, 267, 301, 279, 409, 3, 293, 79, 81, 123, 413, 263, 393, 115, 149, 319, 321, 491, 13, 247, 377, 227, 5, 47, 49, 347, 125, 231, 361, 339, 373, 287, 289, 203, 237, 215, 345, 451, 229, 15, 17, 59, 349, 199, 329, 51, 85, 255, 257, 427, 461, 183, 313, 163, 453, 495, 497, 283, 61, 167, 297, 275, 309, 223, 225, 139, 173, 151, 281, 387, 165, 463, 465, 507, 285, 135, 265, 499, 21, 191, 193, 363, 397, 119, 249, 99, 389, 431, 433, 219, 509, 103, 233, 211, 245, 159, 161, 75, 109, 87, 217, 323, 101, 399, 401, 443, 221, 71, 201, 435, 469, 127, 129, 299, 333, 55, 185, 35, 325, 367, 369, 155, 445, 39, 169, 147, 181, 95, 97, 11, 45, 23, 153, 259, 37, 335, 337, 379, 157, 7, 137, 371, 405, 63, 65, 235, 269, 503, 121, 483, 261, 303, 305, 91, 381, 487, 105, 83, 117, 31, 33, 459, 493, 471, 89, 195, 485, 271, 273, 315, 93, 455, 73, 307, 341, 511};
  static uint16_t N512[512] = {1, 683, 205, 439, 569, 931, 709, 751, 241, 539, 829, 935, 41, 531, 565, 991, 993, 907, 941, 919, 25, 643, 933, 719, 209, 763, 541, 391, 521, 243, 789, 959, 961, 107, 653, 375, 505, 355, 133, 687, 177, 987, 253, 871, 1001, 979, 1013, 927, 929, 331, 365, 855, 985, 67, 357, 655, 145, 187, 989, 327, 457, 691, 213, 895, 897, 555, 77, 311, 441, 803, 581, 623, 113, 411, 701, 807, 937, 403, 437, 863, 865, 779, 813, 791, 921, 515, 805, 591, 81, 635, 413, 263, 393, 115, 661, 831, 833, 1003, 525, 247, 377, 227, 5, 559, 49, 859, 125, 743, 873, 851, 885, 799, 801, 203, 237, 727, 857, 963, 229, 527, 17, 59, 861, 199, 329, 563, 85, 767, 769, 427, 973, 183, 313, 675, 453, 495, 1009, 283, 573, 679, 809, 275, 309, 735, 737, 651, 685, 663, 793, 387, 677, 463, 977, 507, 285, 135, 265, 1011, 533, 703, 705, 875, 397, 119, 249, 99, 901, 431, 945, 731, 1021, 615, 745, 723, 757, 671, 673, 75, 109, 599, 729, 835, 101, 399, 913, 955, 733, 71, 201, 435, 981, 639, 641, 299, 845, 55, 185, 547, 325, 367, 881, 155, 445, 551, 681, 147, 181, 607, 609, 523, 557, 535, 665, 259, 549, 335, 849, 379, 157, 7, 137, 883, 405, 575, 577, 747, 269, 1015, 121, 995, 773, 303, 817, 603, 893, 487, 617, 595, 629, 543, 545, 971, 1005, 471, 601, 707, 997, 271, 785, 827, 605, 967, 73, 307, 853, 511, 513, 171, 717, 951, 57, 419, 197, 239, 753, 27, 317, 423, 553, 19, 53, 479, 481, 395, 429, 407, 537, 131, 421, 207, 721, 251, 29, 903, 9, 755, 277, 447, 449, 619, 141, 887, 1017, 867, 645, 175, 689, 475, 765, 359, 489, 467, 501, 415, 417, 843, 877, 343, 473, 579, 869, 143, 657, 699, 477, 839, 969, 179, 725, 383, 385, 43, 589, 823, 953, 291, 69, 111, 625, 923, 189, 295, 425, 915, 949, 351, 353, 267, 301, 279, 409, 3, 293, 79, 593, 123, 925, 775, 905, 627, 149, 319, 321, 491, 13, 759, 889, 739, 517, 47, 561, 347, 637, 231, 361, 339, 373, 287, 289, 715, 749, 215, 345, 451, 741, 15, 529, 571, 349, 711, 841, 51, 597, 255, 257, 939, 461, 695, 825, 163, 965, 1007, 497, 795, 61, 167, 297, 787, 821, 223, 225, 139, 173, 151, 281, 899, 165, 975, 465, 1019, 797, 647, 777, 499, 21, 191, 193, 363, 909, 631, 761, 611, 389, 943, 433, 219, 509, 103, 233, 211, 245, 159, 161, 587, 621, 87, 217, 323, 613, 911, 401, 443, 221, 583, 713, 947, 469, 127, 129, 811, 333, 567, 697, 35, 837, 879, 369, 667, 957, 39, 169, 659, 693, 95, 97, 11, 45, 23, 153, 771, 37, 847, 337, 891, 669, 519, 649, 371, 917, 63, 65, 235, 781, 503, 633, 483, 261, 815, 305, 91, 381, 999, 105, 83, 117, 31, 33, 459, 493, 983, 89, 195, 485, 783, 273, 315, 93, 455, 585, 819, 341, 1023};
  static uint16_t N1024[1024] = {1, 683, 1229, 1463, 1593, 931, 1733, 1775, 241, 539, 1853, 1959, 1065, 531, 565, 991, 993, 1931, 941, 1943, 1049, 1667, 1957, 1743, 209, 763, 541, 1415, 1545, 243, 1813, 1983, 1985, 1131, 653, 375, 505, 355, 133, 1711, 177, 987, 1277, 871, 2025, 2003, 1013, 927, 929, 331, 365, 855, 2009, 1091, 357, 1679, 145, 1211, 2013, 327, 457, 1715, 213, 1919, 1921, 1579, 77, 1335, 1465, 1827, 581, 1647, 113, 1435, 701, 1831, 937, 1427, 1461, 863, 865, 779, 1837, 1815, 921, 515, 805, 1615, 81, 1659, 1437, 1287, 1417, 1139, 661, 1855, 1857, 2027, 1549, 247, 377, 1251, 1029, 1583, 49, 1883, 125, 743, 1897, 851, 1909, 799, 801, 1227, 1261, 727, 1881, 1987, 1253, 1551, 17, 59, 861, 199, 329, 563, 1109, 1791, 1793, 427, 973, 1207, 1337, 675, 1477, 1519, 2033, 283, 1597, 1703, 809, 275, 309, 735, 737, 1675, 685, 1687, 793, 1411, 1701, 1487, 2001, 507, 285, 1159, 1289, 2035, 1557, 1727, 1729, 875, 397, 119, 249, 99, 1925, 1455, 1969, 731, 1021, 615, 1769, 1747, 757, 671, 673, 75, 109, 599, 1753, 835, 101, 1423, 1937, 955, 1757, 71, 201, 1459, 2005, 1663, 1665, 1323, 1869, 1079, 1209, 1571, 325, 1391, 1905, 1179, 445, 1575, 681, 1171, 1205, 607, 609, 523, 1581, 1559, 665, 259, 549, 1359, 1873, 1403, 1181, 1031, 1161, 883, 405, 1599, 1601, 1771, 1293, 2039, 121, 995, 773, 1327, 1841, 1627, 1917, 487, 1641, 595, 1653, 543, 545, 971, 1005, 471, 1625, 1731, 997, 1295, 1809, 1851, 605, 1991, 73, 307, 853, 1535, 1537, 171, 717, 951, 1081, 419, 1221, 1263, 1777, 27, 1341, 1447, 553, 19, 53, 479, 481, 1419, 429, 1431, 537, 1155, 1445, 1231, 1745, 251, 29, 903, 1033, 1779, 1301, 1471, 1473, 619, 141, 1911, 2041, 1891, 1669, 1199, 1713, 475, 765, 359, 1513, 1491, 501, 415, 417, 1867, 1901, 343, 1497, 579, 1893, 1167, 1681, 699, 1501, 1863, 1993, 1203, 1749, 1407, 1409, 1067, 1613, 823, 953, 1315, 69, 1135, 1649, 923, 189, 1319, 425, 915, 949, 351, 353, 267, 1325, 1303, 409, 3, 293, 1103, 1617, 1147, 925, 775, 905, 627, 149, 1343, 1345, 1515, 1037, 1783, 1913, 739, 517, 1071, 1585, 1371, 1661, 231, 1385, 339, 1397, 287, 289, 715, 749, 215, 1369, 1475, 741, 1039, 1553, 1595, 349, 1735, 1865, 51, 597, 1279, 1281, 1963, 461, 695, 825, 163, 965, 1007, 1521, 1819, 1085, 1191, 297, 1811, 1845, 223, 225, 1163, 173, 1175, 281, 899, 1189, 975, 1489, 2043, 1821, 647, 777, 1523, 1045, 1215, 1217, 363, 1933, 1655, 1785, 1635, 1413, 943, 1457, 219, 509, 103, 1257, 1235, 245, 159, 161, 1611, 1645, 87, 1241, 323, 1637, 911, 1425, 443, 1245, 1607, 1737, 947, 1493, 1151, 1153, 811, 1357, 567, 697, 1059, 1861, 879, 1393, 667, 1981, 1063, 169, 659, 693, 95, 97, 11, 1069, 1047, 153, 1795, 37, 847, 1361, 891, 669, 519, 649, 371, 1941, 1087, 1089, 1259, 781, 1527, 1657, 483, 261, 815, 1329, 1115, 1405, 2023, 1129, 83, 1141, 31, 33, 459, 493, 2007, 1113, 1219, 485, 783, 1297, 1339, 93, 1479, 1609, 1843, 341, 1023, 1025, 1707, 205, 439, 569, 1955, 709, 751, 1265, 1563, 829, 935, 41, 1555, 1589, 2015, 2017, 907, 1965, 919, 25, 643, 933, 719, 1233, 1787, 1565, 391, 521, 1267, 789, 959, 961, 107, 1677, 1399, 1529, 1379, 1157, 687, 1201, 2011, 253, 1895, 1001, 979, 2037, 1951, 1953, 1355, 1389, 1879, 985, 67, 1381, 655, 1169, 187, 989, 1351, 1481, 691, 1237, 895, 897, 555, 1101, 311, 441, 803, 1605, 623, 1137, 411, 1725, 807, 1961, 403, 437, 1887, 1889, 1803, 813, 791, 1945, 1539, 1829, 591, 1105, 635, 413, 263, 393, 115, 1685, 831, 833, 1003, 525, 1271, 1401, 227, 5, 559, 1073, 859, 1149, 1767, 873, 1875, 885, 1823, 1825, 203, 237, 1751, 857, 963, 229, 527, 1041, 1083, 1885, 1223, 1353, 1587, 85, 767, 769, 1451, 1997, 183, 313, 1699, 453, 495, 1009, 1307, 573, 679, 1833, 1299, 1333, 1759, 1761, 651, 1709, 663, 1817, 387, 677, 463, 977, 1531, 1309, 135, 265, 1011, 533, 703, 705, 1899, 1421, 1143, 1273, 1123, 901, 431, 945, 1755, 2045, 1639, 745, 723, 1781, 1695, 1697, 1099, 1133, 1623, 729, 1859, 1125, 399, 913, 1979, 733, 1095, 1225, 435, 981, 639, 641, 299, 845, 55, 185, 547, 1349, 367, 881, 155, 1469, 551, 1705, 147, 181, 1631, 1633, 1547, 557, 535, 1689, 1283, 1573, 335, 849, 379, 157, 7, 137, 1907, 1429, 575, 577, 747, 269, 1015, 1145, 2019, 1797, 303, 817, 603, 893, 1511, 617, 1619, 629, 1567, 1569, 1995, 2029, 1495, 601, 707, 2021, 271, 785, 827, 1629, 967, 1097, 1331, 1877, 511, 513, 1195, 1741, 1975, 57, 1443, 197, 239, 753, 1051, 317, 423, 1577, 1043, 1077, 1503, 1505, 395, 1453, 407, 1561, 131, 421, 207, 721, 1275, 1053, 1927, 9, 755, 277, 447, 449, 1643, 1165, 887, 1017, 867, 645, 175, 689, 1499, 1789, 1383, 489, 467, 1525, 1439, 1441, 843, 877, 1367, 473, 1603, 869, 143, 657, 1723, 477, 839, 969, 179, 725, 383, 385, 43, 589, 1847, 1977, 291, 1093, 111, 625, 1947, 1213, 295, 1449, 1939, 1973, 1375, 1377, 1291, 301, 279, 1433, 1027, 1317, 79, 593, 123, 1949, 1799, 1929, 1651, 1173, 319, 321, 491, 13, 759, 889, 1763, 1541, 47, 561, 347, 637, 1255, 361, 1363, 373, 1311, 1313, 1739, 1773, 1239, 345, 451, 1765, 15, 529, 571, 1373, 711, 841, 1075, 1621, 255, 257, 939, 1485, 1719, 1849, 1187, 1989, 2031, 497, 795, 61, 167, 1321, 787, 821, 1247, 1249, 139, 1197, 151, 1305, 1923, 165, 1999, 465, 1019, 797, 1671, 1801, 499, 21, 191, 193, 1387, 909, 631, 761, 611, 389, 1967, 433, 1243, 1533, 1127, 233, 211, 1269, 1183, 1185, 587, 621, 1111, 217, 1347, 613, 1935, 401, 1467, 221, 583, 713, 1971, 469, 127, 129, 1835, 333, 1591, 1721, 35, 837, 1903, 369, 1691, 957, 39, 1193, 1683, 1717, 1119, 1121, 1035, 45, 23, 1177, 771, 1061, 1871, 337, 1915, 1693, 1543, 1673, 1395, 917, 63, 65, 235, 1805, 503, 633, 1507, 1285, 1839, 305, 91, 381, 999, 105, 1107, 117, 1055, 1057, 1483, 1517, 983, 89, 195, 1509, 1807, 273, 315, 1117, 455, 585, 819, 1365, 2047}; 
  static uint16_t N2048[2048] = {1, 2731, 3277, 3511, 3641, 2979, 3781, 3823, 241, 2587, 3901, 1959, 3113, 2579, 565, 3039, 993, 3979, 2989, 3991, 3097, 3715, 4005, 1743, 2257, 2811, 541, 1415, 3593, 2291, 1813, 4031, 4033, 3179, 653, 2423, 505, 2403, 2181, 1711, 2225, 987, 3325, 2919, 2025, 4051, 1013, 2975, 929, 331, 365, 2903, 4057, 3139, 2405, 3727, 145, 1211, 4061, 2375, 2505, 3763, 2261, 3967, 3969, 3627, 2125, 1335, 1465, 1827, 581, 3695, 113, 3483, 2749, 3879, 937, 1427, 1461, 2911, 865, 779, 1837, 1815, 921, 2563, 805, 1615, 2129, 3707, 3485, 3335, 1417, 1139, 2709, 3903, 3905, 4075, 3597, 247, 2425, 1251, 3077, 1583, 2097, 1883, 2173, 743, 3945, 2899, 1909, 2847, 801, 1227, 3309, 727, 1881, 1987, 3301, 3599, 17, 2107, 2909, 199, 329, 2611, 3157, 3839, 3841, 427, 973, 3255, 3385, 675, 1477, 3567, 4081, 283, 1597, 1703, 2857, 275, 2357, 2783, 737, 1675, 685, 3735, 2841, 1411, 1701, 1487, 2001, 507, 2333, 1159, 3337, 4083, 3605, 3775, 3777, 875, 2445, 2167, 249, 99, 3973, 1455, 1969, 2779, 1021, 2663, 1769, 1747, 2805, 2719, 673, 2123, 2157, 2647, 3801, 835, 101, 3471, 3985, 3003, 1757, 2119, 2249, 1459, 4053, 3711, 3713, 1323, 3917, 1079, 1209, 3619, 2373, 3439, 3953, 1179, 445, 3623, 681, 3219, 3253, 2655, 609, 2571, 3629, 1559, 665, 259, 2597, 1359, 1873, 1403, 1181, 3079, 1161, 2931, 405, 3647, 3649, 1771, 1293, 4087, 2169, 3043, 773, 1327, 1841, 3675, 3965, 487, 3689, 595, 3701, 2591, 545, 3019, 1005, 471, 1625, 3779, 997, 3343, 3857, 3899, 605, 4039, 73, 307, 853, 3583, 3585, 2219, 2765, 2999, 3129, 2467, 3269, 3311, 3825, 2075, 3389, 1447, 2601, 2067, 53, 2527, 481, 3467, 2477, 3479, 2585, 3203, 3493, 1231, 1745, 2299, 29, 903, 3081, 1779, 1301, 3519, 3521, 2667, 141, 1911, 4089, 1891, 1669, 1199, 1713, 475, 2813, 2407, 1513, 3539, 501, 2463, 417, 3915, 3949, 2391, 3545, 2627, 1893, 3215, 3729, 699, 3549, 1863, 1993, 3251, 1749, 3455, 3457, 3115, 1613, 823, 953, 1315, 69, 3183, 3697, 2971, 2237, 3367, 425, 915, 949, 2399, 353, 267, 1325, 1303, 409, 2051, 293, 1103, 1617, 3195, 2973, 2823, 905, 627, 2197, 3391, 3393, 3563, 3085, 3831, 1913, 739, 2565, 1071, 1585, 1371, 1661, 231, 3433, 2387, 1397, 2335, 289, 715, 2797, 215, 1369, 1475, 2789, 3087, 3601, 1595, 2397, 3783, 3913, 2099, 2645, 3327, 3329, 4011, 461, 2743, 2873, 163, 965, 3055, 3569, 3867, 1085, 1191, 2345, 3859, 1845, 2271, 225, 1163, 173, 3223, 2329, 899, 1189, 975, 1489, 4091, 1821, 647, 2825, 3571, 3093, 3263, 3265, 363, 1933, 1655, 3833, 3683, 3461, 943, 1457, 2267, 509, 2151, 1257, 1235, 2293, 2207, 161, 1611, 1645, 2135, 3289, 323, 3685, 2959, 3473, 2491, 1245, 1607, 1737, 947, 3541, 3199, 3201, 811, 3405, 567, 697, 3107, 1861, 2927, 3441, 667, 4029, 3111, 169, 2707, 2741, 2143, 97, 2059, 3117, 1047, 153, 3843, 2085, 847, 1361, 891, 669, 2567, 649, 2419, 3989, 3135, 3137, 1259, 781, 3575, 1657, 2531, 261, 815, 1329, 3163, 3453, 4071, 3177, 83, 3189, 2079, 33, 2507, 493, 4055, 1113, 3267, 485, 2831, 3345, 3387, 93, 3527, 3657, 3891, 341, 3071, 3073, 1707, 2253, 2487, 2617, 1955, 2757, 2799, 3313, 1563, 2877, 935, 2089, 1555, 3637, 2015, 4065, 2955, 1965, 2967, 2073, 2691, 2981, 719, 1233, 1787, 3613, 391, 2569, 1267, 789, 3007, 3009, 2155, 3725, 1399, 3577, 1379, 1157, 687, 1201, 4059, 2301, 1895, 1001, 3027, 4085, 1951, 4001, 3403, 3437, 1879, 3033, 2115, 1381, 2703, 3217, 187, 3037, 1351, 1481, 2739, 1237, 2943, 2945, 2603, 1101, 311, 441, 803, 3653, 2671, 3185, 2459, 1725, 2855, 4009, 403, 437, 1887, 3937, 3851, 813, 791, 3993, 1539, 3877, 591, 1105, 2683, 2461, 2311, 393, 115, 1685, 2879, 2881, 3051, 2573, 3319, 1401, 227, 2053, 559, 1073, 859, 1149, 3815, 2921, 1875, 885, 1823, 3873, 203, 2285, 3799, 857, 963, 2277, 2575, 3089, 1083, 1885, 3271, 3401, 1587, 2133, 2815, 2817, 3499, 4045, 2231, 2361, 3747, 453, 2543, 3057, 3355, 573, 679, 1833, 3347, 1333, 1759, 3809, 651, 3757, 2711, 1817, 387, 677, 463, 977, 3579, 1309, 135, 2313, 3059, 2581, 2751, 2753, 3947, 1421, 1143, 3321, 3171, 2949, 431, 945, 1755, 4093, 1639, 745, 723, 1781, 1695, 3745, 1099, 1133, 1623, 2777, 3907, 3173, 2447, 2961, 1979, 733, 1095, 1225, 435, 3029, 2687, 2689, 299, 2893, 55, 185, 2595, 1349, 2415, 2929, 155, 3517, 2599, 3753, 2195, 2229, 1631, 3681, 1547, 2605, 535, 3737, 3331, 1573, 335, 849, 379, 157, 2055, 137, 1907, 3477, 2623, 2625, 747, 269, 3063, 1145, 2019, 3845, 303, 817, 2651, 2941, 3559, 2665, 3667, 2677, 1567, 3617, 1995, 4077, 3543, 601, 2755, 4069, 2319, 2833, 2875, 3677, 3015, 3145, 3379, 3925, 2559, 2561, 1195, 1741, 1975, 2105, 1443, 2245, 2287, 2801, 1051, 2365, 423, 1577, 1043, 3125, 1503, 3553, 2443, 1453, 2455, 1561, 2179, 2469, 207, 721, 1275, 3101, 3975, 2057, 755, 277, 2495, 2497, 1643, 3213, 887, 3065, 867, 645, 175, 689, 3547, 1789, 1383, 489, 2515, 3573, 1439, 3489, 2891, 2925, 1367, 2521, 1603, 869, 2191, 2705, 3771, 2525, 839, 969, 2227, 725, 2431, 2433, 2091, 589, 3895, 4025, 291, 3141, 2159, 2673, 1947, 1213, 2343, 3497, 3987, 4021, 1375, 3425, 3339, 301, 279, 3481, 1027, 3365, 79, 593, 2171, 1949, 1799, 3977, 3699, 1173, 2367, 2369, 2539, 2061, 2807, 889, 3811, 1541, 47, 561, 347, 637, 3303, 2409, 1363, 373, 1311, 3361, 3787, 1773, 3287, 345, 451, 1765, 2063, 2577, 571, 1373, 2759, 2889, 1075, 1621, 2303, 2305, 2987, 3533, 1719, 1849, 3235, 4037, 2031, 2545, 2843, 61, 167, 1321, 2835, 821, 1247, 3297, 139, 3245, 2199, 1305, 3971, 165, 4047, 465, 3067, 797, 3719, 1801, 2547, 2069, 2239, 2241, 3435, 909, 631, 2809, 2659, 2437, 4015, 433, 1243, 3581, 1127, 233, 211, 1269, 1183, 3233, 587, 621, 1111, 2265, 3395, 2661, 1935, 2449, 1467, 221, 583, 713, 4019, 2517, 2175, 2177, 3883, 2381, 3639, 3769, 2083, 837, 1903, 2417, 3739, 3005, 2087, 3241, 1683, 1717, 1119, 3169, 1035, 2093, 23, 3225, 2819, 1061, 3919, 337, 3963, 3741, 1543, 3721, 1395, 2965, 2111, 2113, 235, 3853, 2551, 633, 1507, 3333, 3887, 305, 2139, 2429, 3047, 2153, 3155, 2165, 1055, 3105, 1483, 3565, 3031, 89, 2243, 3557, 1807, 2321, 2363, 3165, 2503, 2633, 2867, 3413, 2047, 2049, 683, 1229, 1463, 1593, 931, 1733, 1775, 2289, 539, 1853, 4007, 1065, 531, 2613, 991, 3041, 1931, 941, 1943, 1049, 1667, 1957, 3791, 209, 763, 2589, 3463, 1545, 243, 3861, 1983, 1985, 1131, 2701, 375, 2553, 355, 133, 3759, 177, 3035, 1277, 871, 4073, 2003, 3061, 927, 2977, 2379, 2413, 855, 2009, 1091, 357, 1679, 2193, 3259, 2013, 327, 457, 1715, 213, 1919, 1921, 1579, 77, 3383, 3513, 3875, 2629, 1647, 2161, 1435, 701, 1831, 2985, 3475, 3509, 863, 2913, 2827, 3885, 3863, 2969, 515, 2853, 3663, 81, 1659, 1437, 1287, 3465, 3187, 661, 1855, 1857, 2027, 1549, 2295, 377, 3299, 1029, 3631, 49, 3931, 125, 2791, 1897, 851, 3957, 799, 2849, 3275, 1261, 2775, 3929, 4035, 1253, 1551, 2065, 59, 861, 2247, 2377, 563, 1109, 1791, 1793, 2475, 3021, 1207, 1337, 2723, 3525, 1519, 2033, 2331, 3645, 3751, 809, 2323, 309, 735, 2785, 3723, 2733, 1687, 793, 3459, 3749, 3535, 4049, 2555, 285, 3207, 1289, 2035, 1557, 1727, 1729, 2923, 397, 119, 2297, 2147, 1925, 3503, 4017, 731, 3069, 615, 3817, 3795, 757, 671, 2721, 75, 109, 599, 1753, 2883, 2149, 1423, 1937, 955, 3805, 71, 201, 3507, 2005, 1663, 1665, 3371, 1869, 3127, 3257, 1571, 325, 1391, 1905, 3227, 2493, 1575, 2729, 1171, 1205, 607, 2657, 523, 1581, 3607, 2713, 2307, 549, 3407, 3921, 3451, 3229, 1031, 3209, 883, 2453, 1599, 1601, 3819, 3341, 2039, 121, 995, 2821, 3375, 3889, 1627, 1917, 2535, 1641, 2643, 1653, 543, 2593, 971, 3053, 2519, 3673, 1731, 3045, 1295, 1809, 1851, 2653, 1991, 2121, 2355, 2901, 1535, 1537, 171, 717, 951, 1081, 419, 1221, 1263, 1777, 27, 1341, 3495, 553, 19, 2101, 479, 2529, 1419, 429, 1431, 537, 1155, 1445, 3279, 3793, 251, 2077, 2951, 1033, 3827, 3349, 1471, 1473, 619, 2189, 3959, 2041, 3939, 3717, 3247, 3761, 2523, 765, 359, 3561, 1491, 2549, 415, 2465, 1867, 1901, 343, 1497, 579, 3941, 1167, 1681, 2747, 1501, 3911, 4041, 1203, 3797, 1407, 1409, 1067, 3661, 2871, 3001, 3363, 2117, 1135, 1649, 923, 189, 1319, 2473, 2963, 2997, 351, 2401, 2315, 3373, 3351, 2457, 3, 2341, 3151, 3665, 1147, 925, 775, 2953, 2675, 149, 1343, 1345, 1515, 1037, 1783, 3961, 2787, 517, 3119, 3633, 3419, 3709, 2279, 1385, 339, 3445, 287, 2337, 2763, 749, 2263, 3417, 3523, 741, 1039, 1553, 3643, 349, 1735, 1865, 51, 597, 1279, 1281, 1963, 2509, 695, 825, 2211, 3013, 1007, 1521, 1819, 3133, 3239, 297, 1811, 3893, 223, 2273, 3211, 2221, 1175, 281, 2947, 3237, 3023, 3537, 2043, 3869, 2695, 777, 1523, 1045, 1215, 1217, 2411, 3981, 3703, 1785, 1635, 1413, 2991, 3505, 219, 2557, 103, 3305, 3283, 245, 159, 2209, 3659, 3693, 87, 1241, 2371, 1637, 911, 1425, 443, 3293, 3655, 3785, 2995, 1493, 1151, 1153, 2859, 1357, 2615, 2745, 1059, 3909, 879, 1393, 2715, 1981, 1063, 2217, 659, 693, 95, 2145, 11, 1069, 3095, 2201, 1795, 37, 2895, 3409, 2939, 2717, 519, 2697, 371, 1941, 1087, 1089, 3307, 2829, 1527, 3705, 483, 2309, 2863, 3377, 1115, 1405, 2023, 1129, 2131, 1141, 31, 2081, 459, 2541, 2007, 3161, 1219, 2533, 783, 1297, 1339, 2141, 1479, 1609, 1843, 2389, 1023, 1025, 3755, 205, 439, 569, 4003, 709, 751, 1265, 3611, 829, 2983, 41, 3603, 1589, 4063, 2017, 907, 4013, 919, 25, 643, 933, 2767, 3281, 3835, 1565, 2439, 521, 3315, 2837, 959, 961, 107, 1677, 3447, 1529, 3427, 3205, 2735, 3249, 2011, 253, 3943, 3049, 979, 2037, 3999, 1953, 1355, 1389, 3927, 985, 67, 3429, 655, 1169, 2235, 989, 3399, 3529, 691, 3285, 895, 897, 555, 3149, 2359, 2489, 2851, 1605, 623, 1137, 411, 3773, 807, 1961, 2451, 2485, 3935, 1889, 1803, 2861, 2839, 1945, 3587, 1829, 2639, 3153, 635, 413, 263, 2441, 2163, 3733, 831, 833, 1003, 525, 1271, 3449, 2275, 5, 2607, 3121, 2907, 3197, 1767, 873, 3923, 2933, 3871, 1825, 2251, 237, 1751, 2905, 3011, 229, 527, 1041, 3131, 3933, 1223, 1353, 3635, 85, 767, 769, 1451, 1997, 183, 313, 1699, 2501, 495, 1009, 1307, 2621, 2727, 3881, 1299, 3381, 3807, 1761, 2699, 1709, 663, 3865, 2435, 2725, 2511, 3025, 1531, 3357, 2183, 265, 1011, 533, 703, 705, 1899, 3469, 3191, 1273, 1123, 901, 2479, 2993, 3803, 2045, 3687, 2793, 2771, 3829, 3743, 1697, 3147, 3181, 3671, 729, 1859, 1125, 399, 913, 4027, 2781, 3143, 3273, 2483, 981, 639, 641, 2347, 845, 2103, 2233, 547, 3397, 367, 881, 2203, 1469, 551, 1705, 147, 181, 3679, 1633, 3595, 557, 2583, 1689, 1283, 3621, 2383, 2897, 2427, 2205, 7, 2185, 3955, 1429, 575, 577, 2795, 2317, 1015, 3193, 4067, 1797, 2351, 2865, 603, 893, 1511, 617, 1619, 629, 3615, 1569, 4043, 2029, 1495, 2649, 707, 2021, 271, 785, 827, 1629, 967, 1097, 1331, 1877, 511, 513, 3243, 3789, 4023, 57, 3491, 197, 239, 753, 3099, 317, 2471, 3625, 3091, 1077, 3551, 1505, 395, 3501, 407, 3609, 131, 421, 2255, 2769, 3323, 1053, 1927, 9, 2803, 2325, 447, 449, 3691, 1165, 2935, 1017, 2915, 2693, 2223, 2737, 1499, 3837, 3431, 2537, 467, 1525, 3487, 1441, 843, 877, 3415, 473, 3651, 2917, 143, 657, 1723, 477, 2887, 3017, 179, 2773, 383, 385, 43, 2637, 1847, 1977, 2339, 1093, 111, 625, 3995, 3261, 295, 1449, 1939, 1973, 3423, 1377, 1291, 2349, 2327, 1433, 3075, 1317, 2127, 2641, 123, 3997, 3847, 1929, 1651, 3221, 319, 321, 491, 13, 759, 2937, 1763, 3589, 2095, 2609, 2395, 2685, 1255, 361, 3411, 2421, 3359, 1313, 1739, 3821, 1239, 2393, 2499, 3813, 15, 529, 2619, 3421, 711, 841, 3123, 3669, 255, 257, 939, 1485, 3767, 3897, 1187, 1989, 4079, 497, 795, 2109, 2215, 3369, 787, 2869, 3295, 1249, 2187, 1197, 151, 3353, 1923, 2213, 1999, 2513, 1019, 2845, 1671, 3849, 499, 21, 191, 193, 1387, 2957, 2679, 761, 611, 389, 1967, 2481, 3291, 1533, 3175, 2281, 2259, 3317, 3231, 1185, 2635, 2669, 3159, 217, 1347, 613, 3983, 401, 3515, 2269, 2631, 2761, 1971, 469, 127, 129, 1835, 333, 1591, 1721, 35, 2885, 3951, 369, 1691, 957, 39, 1193, 3731, 3765, 3167, 1121, 3083, 45, 2071, 1177, 771, 3109, 1871, 2385, 1915, 1693, 3591, 1673, 3443, 917, 63, 65, 2283, 1805, 503, 2681, 3555, 1285, 1839, 2353, 91, 381, 999, 105, 1107, 117, 3103, 1057, 3531, 1517, 983, 2137, 195, 1509, 3855, 273, 315, 1117, 455, 585, 819, 1365, 4095};
  static uint16_t * N4096 = NULL;
  static uint16_t * N_precalc[9] = {N256, N512, N1024, NULL, N2048, NULL, NULL, NULL, NULL};
  if(!N4096) {
    N4096 = pre_compute_odd_inverse_mod(2*4096);
    N_precalc[8] = N4096;
  }
  return N_precalc[N>>9][(x-1)>>1];
}

// generic inverse mod
// From hexl https://github.com/intel/hexl
// Apache-2.0 license - Copyright 2020-2021 Intel Corporation
uint64_t inverse_mod(uint64_t input, uint64_t modulus) {
  uint64_t a = input % modulus;
  assert(a != 0);
  if (modulus == 1) {
    return 0;
  }

  int64_t m0 = (int64_t)(modulus);
  int64_t y = 0;
  int64_t x = 1;
  while (a > 1) {
    // q is quotient
    int64_t q = (int64_t)(a / modulus);

    int64_t t = (int64_t)(modulus);
    modulus = a % modulus;
    a = (uint64_t)(t);

    // Update y and x
    t = y;
    y = x - q * y;
    x = t;
  }

  // Make x positive
  if (x < 0) x += m0;

  return (uint64_t)(x);
}


uint64_t * pre_compute_inverse_mod(uint64_t modulus){
  uint64_t * res = (uint64_t *) safe_aligned_malloc(sizeof(uint64_t)*modulus-1);
  for (size_t i = 0; i < modulus; i++){
    res[i] = inverse_mod(i, modulus); 
  }
  return res;
}