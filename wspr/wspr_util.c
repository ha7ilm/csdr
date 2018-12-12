/*
 This file is part of program wsprd, a detector/demodulator/decoder
 for the Weak Signal Propagation Reporter (WSPR) mode.
 
 File name: wsprd_util.c
 
 Copyright 2001-2015, Joe Taylor, K1JT
 
 Most of the code is based on work by Steven Franke, K9AN, which
 in turn was based on earlier work by K1JT.
 
 Copyright 2014-2015, Steven Franke, K9AN
 
 License: GNU GPL v3
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "wspr.h"
#include "nhash.h"

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <strings.h>
#include <sys/time.h>

// call_28b bits 27-0 = first 28 bits of *d (big-endian)
// ......d0 ......d1 ......d2 ......d3
// 22222222 11111111 11000000 0000....
// 76543210 98765432 10987654 3210....

// grid_pwr_22b bits 21-0 = 22 subsequent bits
// ......d3 ......d4 ......d5 ......d6
// ....2211 11111111 00000000 00......
// ....1098 76543210 98765432 10......
//     gggg gggggggg gggppppp pp		g: 15-bit grid, p: 7-bit pwr

void unpack50(u1_t *d, u4_t *call_28b, u4_t *grid_pwr_22b, u4_t *grid_15b, u4_t *pwr_7b)
{
    *call_28b = (d[0]<<20) | (d[1]<<12) | (d[2]<<4) | (d[3]>>4);
    *grid_pwr_22b = ((d[3]&0xf)<<18) | (d[4]<<10) | (d[5]<<2) | (d[6]>>6);
	*grid_15b = *grid_pwr_22b >> 7;
	*pwr_7b = *grid_pwr_22b & 0x7f;
}

static const char *c = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ ";		// 37 characters

// valid call: AANLLL, A = alpha-numeric, N = numeric, L = letter or space
int unpackcall(u4_t call_28b, char *call)
{
    int i;
    char tmp[7];
    
    memset(call, 0, LEN_CALL);
    
    if (call_28b < 262177560) {		// = 37*36*10*27*27*27 = 0xfa08318 (28-bits)
        tmp[5] = c[call_28b%27+10];	// letter, space
        call_28b /= 27;
        tmp[4] = c[call_28b%27+10];	// letter, space
        call_28b /= 27;
        tmp[3] = c[call_28b%27+10];	// letter, space
        call_28b /= 27;
        tmp[2] = c[call_28b%10];	// number
        call_28b /= 10;
        tmp[1] = c[call_28b%36];	// letter, number
        call_28b /= 36;
        tmp[0] = c[call_28b];		// letter, number, space
        tmp[6] = '\0';

		// remove leading whitespace
        for (i=0; i<5; i++) {
            if (tmp[i] != ' ')
                break;
        }
        sprintf(call, "%-6s", &tmp[i]);
        
		// remove trailing whitespace
        for (i=0; i<6; i++) {
            if (call[i] == ' ') {
                call[i] = '\0';
            }
        }
    } else {
    	return 0;
    }
    return 1;
}

int unpackgrid(u4_t grid_15b, char *grid)
{
    int dlat, dlong;
    
    memset(grid, 0, LEN_GRID);

    if (grid_15b < 32400) {		// = 0x7e90 (15-bits)
        dlat = (grid_15b%180)-90;
        dlong = (grid_15b/180)*2 - 180 + 2;
        if (dlong < -180)
            dlong = dlong+360;
        if (dlong > 180)
            dlong = dlong+360;
        int nlong = 60.0*(180.0-dlong)/5.0;
        int n1 = nlong/240;
        int n2 = (nlong - 240*n1)/24;
        //int n3 = nlong -40*n1 - 24*n2;
        grid[0] = c[10+n1];
        grid[2] = c[n2];

        int nlat = 60.0*(dlat+90)/2.5;
        n1 = nlat/240;
        n2 = (nlat-240*n1)/24;
        //n3 = nlong - 240*n1 - 24*n2;
        grid[1] = c[10+n1];
        grid[3] = c[n2];
    } else {
        strcpy(grid,"XXXX");
        return 0;
    }
    return 1;
}

int unpackpfx(int32_t nprefix, char *call)
{
    char nc, pfx[4]={'\0'}, tmpcall[7];
    int i;
    int32_t n;
    
    strcpy(tmpcall,call);
    if( nprefix < 60000 ) {		// < 0xea60
        // add a prefix of 1 to 3 characters
        n=nprefix;
        for (i=2; i>=0; i--) {
            nc = n%37;		// letter, number, space
            if( (nc >= 0) & (nc <= 9) ) {
                pfx[i] = nc + '0';
            }
            else if( (nc >= 10) & (nc <= 35) ) {
                pfx[i] = nc-10 + 'A';
            }
            else {
                pfx[i] = ' ';
            }
            n = n/37;
        }

        char * p = strrchr(pfx,' ');
        strcpy(call, p ? p + 1 : pfx);
        strncat(call,"/",1);
        strncat(call,tmpcall,strlen(tmpcall));
        
    } else {
        // add a suffix of 1 or 2 characters
        nc=nprefix-60000;
        if( (nc >= 0) & (nc <= 9) ) {
            pfx[0] = nc + '0';
            strcpy(call,tmpcall);
            strncat(call,"/",1);
            strncat(call,pfx,1);
        }
        else if( (nc >= 10) & (nc <= 35) ) {
            pfx[0] = nc-10 + 'A';
            strcpy(call,tmpcall);
            strncat(call,"/",1);
            strncat(call,pfx,1);
        }
        else if( (nc >= 36) & (nc <= 125) ) {
            pfx[0]=(nc-26)/10 + '0';
            pfx[1]=(nc-26)%10 + '0';
            strcpy(call,tmpcall);
            strncat(call,"/",1);
            strncat(call,pfx,2);
        }
        else {
            return 0;
        }
    }
    return 1;
}

void deinterleave(unsigned char *sym)
{
    unsigned char tmp[NSYM_162];
    unsigned char p, i, j;
    
    p=0;
    i=0;
    while (p<NSYM_162) {
        j=((i * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> 32;
        if (j < NSYM_162 ) {
            tmp[p]=sym[j];
            p=p+1;
        }
        i=i+1;
    }
    for (i=0; i<NSYM_162; i++)
        sym[i]=tmp[i];
}

#define WSPR_HASH_ENTRY_SIZE 16

typedef struct {
	u2_t hash;
	union {
		char call[LEN_CALL];
		char pad[WSPR_HASH_ENTRY_SIZE - sizeof(u2_t)];
	};
} hashtab_t;

static hashtab_t *ht;
static int htsize = 16;

void wspr_hash_init()
{
	assert(sizeof(hashtab_t) == WSPR_HASH_ENTRY_SIZE);
	assert(LEN_CALL <= (WSPR_HASH_ENTRY_SIZE - sizeof(u2_t)));
	int i;
	ht = (hashtab_t *) calloc(htsize, sizeof(hashtab_t));
	assert(ht != NULL);
}

void hash_update(char *call)
{
	int i;
	u2_t hash = nhash(call, strlen(call), (uint32_t) 146);
	
	for (i=0; i < htsize; i++) {
		if (ht[i].call[0] == '\0') {
			//wprintf("W-HASH %d 0x%04x upd new %s\n", i, hash, call);
			ht[i].hash = hash;
			strcpy(ht[i].call, call);
			break;
		}
		if (ht[i].hash == hash) {
			if (strcmp(ht[i].call, call) == 0) {
				//wprintf("W-HASH %d 0x%04x upd hit %s\n", i, hash, call);
			} else {
				//wprintf("W-HASH %d 0x%04x upd COLLISION %s %s\n", i, hash, ht[i].call, call);
				strcpy(ht[i].call, call);
			}
			break;
		}
	}
	
	if (i == htsize) {
		//wprintf("W-HASH expand %d -> %d\n", htsize, htsize*2);
		htsize *= 2;
		ht = (hashtab_t *) realloc(ht, sizeof(hashtab_t) * htsize);
		memset(ht + htsize/2, 0, sizeof(hashtab_t) * htsize/2);
		//wprintf("W-HASH %d 0x%04x exp new %s\n", htsize/2, hash, call);
		ht[htsize/2].hash = hash;
		strcpy(ht[htsize/2].call, call);
	}
}

static char *hash_lookup(int hash)
{
	int i;
	
	for (i=0; i < htsize; i++) {
		if (ht[i].call[0] == '\0')
			break;
		if (ht[i].hash == hash) {
			//wprintf("W-HASH %d 0x%04x lookup %s\n", i, hash, ht[i].call);
			return ht[i].call;
		}
	}
	
	//wprintf("W-HASH 0x%04x lookup FAIL\n", hash);
	return NULL;
}

int unpk_(u1_t *decdata, char *call_loc_pow, char *callsign, char *grid, int *dBm)
{
	int rtn = 0;
    u4_t call_28b, grid_pwr_22b, grid_15b, pwr_7b;
    int n3, ndbm, nadd;
    
    memset(call_loc_pow, 0, LEN_C_L_P);

    unpack50(decdata, &call_28b, &grid_pwr_22b, &grid_15b, &pwr_7b);
    if (!unpackcall(call_28b, callsign)) return -1;
    if (!unpackgrid(grid_15b, grid)) return -2;

	// ntype -64..63, but this is NOT twos complement (just biasing)
    int ntype = pwr_7b - 64;

    /*
     Based on the value of ntype, decide whether this is a Type 1, 2, or 3 message.
     
     * Type 1: 6 digit call, grid, power - ntype is positive and is a member
     of the set {0,3,7,10,13,17,20...60}
     
     * Type 2: extended callsign, power - ntype is positive but not
     a member of the set of allowed powers
     
     * Type 3: hash, 6 digit grid, power - ntype is negative.
     */

    if ((ntype >= 0) && (ntype <= 62)) {
        int nu = ntype%10;
        if (nu == 0 || nu == 3 || nu == 7) {
            *dBm = ndbm = ntype;
        	sprintf(call_loc_pow, "%s %s %2d", callsign, grid, ndbm);
			hash_update(callsign);
            rtn = 1;
        } else {
            nadd = nu;
            if( nu > 3 ) nadd=nu-3;
            if( nu > 7 ) nadd=nu-7;
            // Nggggggg gggggggg
            n3 = grid_15b + 32768*(nadd-1);		// 32768 = 0x8000, (nadd-1) = 0..2
            if (!unpackpfx(n3, callsign)) return -3;
            
            grid[0] = '\0';		// no grid info for TYPE2

            *dBm = ndbm = ntype-nadd;
        	sprintf(call_loc_pow, "%s no_grid %2d", callsign, ndbm);

            int nu = ndbm%10;
            if (nu == 0 || nu == 3 || nu == 7 || nu == 10) { // make sure power is OK
				hash_update(callsign);
            } else return -4;
            rtn = 2;
        }
    } else

    if (ntype < 0) {
        *dBm = ndbm = -(ntype+1);

        memset(grid, 0, LEN_GRID);
        // this because the grid6 L1L2N1N2L3L4 was packed as L2N1N2L3L4L1 to be
        // compatible with the AANLLL required by pack_call() / unpackcall()
        sprintf(grid, "%c%.5s", callsign[5], callsign);

        int nu = ndbm%10;
        if ((nu != 0 && nu != 3 && nu != 7 && nu != 10) ||
            !isalpha(grid[0]) || !isalpha(grid[1]) ||
            !isdigit(grid[2]) || !isdigit(grid[3])) {
               // not testing 4'th and 5'th chars because of this case: <PA0SKT/2> JO33 40
               // grid is only 4 chars even though this is a hashed callsign...
               //         isalpha(grid[4]) && isalpha(grid[5]) ) ) {
        	return -5;
        }
        
        // get hashed callsign
        // wspr_t:callsign shouldn't contain HTML-confusing angle brackets
        char *callp = hash_lookup((grid_pwr_22b-ntype-64) / 128);
        sprintf(callsign, "%s", callp? callp : "...");
        sprintf(call_loc_pow, "<%s> %s %2d", callsign, grid, ndbm);
        
        // I don't know what to do with these... They show up as "A000AA" grids.
        if (ntype == -64) return -6;
        rtn = 3;
    } else {
    	return -7;
    }
    return rtn;
}

int snr_comp(const void *elem1, const void *elem2)
{
	const pk_t *e1 = (const pk_t*) elem1, *e2 = (const pk_t*) elem2;
	int r = (e1->snr0 < e2->snr0)? 1 : ((e1->snr0 > e2->snr0)? -1:0);	// NB: comparison reversed to sort in descending order
	return r;
}

int freq_comp(const void *elem1, const void *elem2)
{
	const pk_t *e1 = (const pk_t*) elem1, *e2 = (const pk_t*) elem2;
	int r = (e1->freq0 < e2->freq0)? -1 : ((e1->freq0 > e2->freq0)? 1:0);
	return r;
}

#define DEG_2_RAD(deg) ((deg) * K_PI / 180.0)

#define latLon_deg_to_rad(loc) \
	loc.lat = DEG_2_RAD(loc.lat); \
	loc.lon = DEG_2_RAD(loc.lon);

// field square subsquare (extended square)
//   A-R    0-9       a-x              0-9
//   #18    #10       #24              #10

#define SQ_LON_DEG		2.0
#define SQ_LAT_DEG		1.0
#define SUBSQ_PER_SQ	24.0
#define SUBSQ_LON_DEG	(SQ_LON_DEG / SUBSQ_PER_SQ)
#define SUBSQ_LAT_DEG	(SQ_LAT_DEG / SUBSQ_PER_SQ)

#define SQ_PER_FLD		10.0
#define	FLD_DEG_LON		(SQ_PER_FLD * SQ_LON_DEG)
#define	FLD_DEG_LAT		(SQ_PER_FLD * SQ_LAT_DEG)

static void grid_to_latLon(char *grid, latLon_t *loc)
{
	double lat, lon;
	char c;
	int slen = strlen(grid);
	
	loc->lat = loc->lon = 999.0;
	if (slen < 4) return;
	
	c = tolower(grid[0]);
	if (c < 'a' || c > 'r') return;
	lon = (c-'a')*20 - 180;

	c = tolower(grid[1]);
	if (c < 'a' || c > 'r') return;
	lat = (c-'a')*10 - 90;

	c = grid[2];
	if (c < '0' || c > '9') return;
	lon += (c-'0') * SQ_LON_DEG;

	c = grid[3];
	if (c < '0' || c > '9') return;
	lat += (c-'0') * SQ_LAT_DEG;

	if (slen != 6) {	// assume center of square (i.e. "....ll")
		lon += SQ_LON_DEG /2.0;
		lat += SQ_LAT_DEG /2.0;
	} else {
		c = tolower(grid[4]);
		if (c < 'a' || c > 'x') return;
		lon += (c-'a') * SUBSQ_LON_DEG;

		c = tolower(grid[5]);
		if (c < 'a' || c > 'x') return;
		lat += (c-'a') * SUBSQ_LAT_DEG;

		lon += SUBSQ_LON_DEG /2.0;	// assume center of sub-square (i.e. "......44")
		lat += SUBSQ_LAT_DEG /2.0;
	}

	loc->lat = lat;
	loc->lon = lon;
	//wprintf("GRID %s%s = (%f, %f)\n", grid, (slen != 6)? "[ll]":"", lat, lon);
}

static const char *field = "ABCDEFGHIJKLMNOPQR";
static const char *square = "0123456789";
static const char *subsquare = "abcdefghijklmnopqrstuvwx";

int latLon_to_grid6(latLon_t *loc, char *grid6)
{
	int i;
	double r, lat, lon;
	
	// longitude
	lon = loc->lon + 180.0;
	if (lon < 0 || lon >= 360.0) return -1;
	i = (int) lon / FLD_DEG_LON;
	grid6[0] = field[i];
	r = lon - (i * FLD_DEG_LON);
	
	i = (int) floor(r / SQ_LON_DEG);
	grid6[2] = square[i];
	r = r - (i * SQ_LON_DEG);
	
	i = (int) floor(r * (SUBSQ_PER_SQ / SQ_LON_DEG));
	grid6[4] = subsquare[i];
	
	// latitude
	lat = loc->lat + 90.0;
	if (lat < 0 || lat >= 180.0) return -1;
	i = (int) lat / FLD_DEG_LAT;
	grid6[1] = field[i];
	r = lat - (i * FLD_DEG_LAT);
	
	i = (int) floor(r / SQ_LAT_DEG);
	grid6[3] = square[i];
	r = r - (i * SQ_LAT_DEG);
	
	i = (int) floor(r * (SUBSQ_PER_SQ / SQ_LAT_DEG));
	grid6[5] = subsquare[i];
	
	return 0;
}

static latLon_t r_loc;

void set_reporter_grid(char *grid)
{
	grid_to_latLon(grid, &r_loc);
	if (r_loc.lat != 999.0)
		latLon_deg_to_rad(r_loc);
}

double grid_to_distance_km(char *grid)
{
	if (r_loc.lat == 999.0 || *grid == '\0')
		return 0;
	
	latLon_t loc;
	grid_to_latLon(grid, &loc);
	latLon_deg_to_rad(loc);
	
	double delta_lat = loc.lat - r_loc.lat;
	delta_lat /= 2.0;
	delta_lat = sin(delta_lat);
	delta_lat *= delta_lat;
	double delta_lon = loc.lon - r_loc.lon;
	delta_lon /= 2.0;
	delta_lon = sin(delta_lon);
	delta_lon *= delta_lon;

	double t = delta_lat + (delta_lon * cos(loc.lat) * cos(r_loc.lat));
	#define EARTH_RADIUS_KM 6371.0
	double km = EARTH_RADIUS_KM * 2.0 * atan2(sqrt(t), sqrt(1.0-t));
	return km;
}
