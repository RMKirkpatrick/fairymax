/***************************************************************************/
/*                               Fairy-Max,                                */
/* Version of the sub-2KB (source) micro-Max Chess program, fused to a     */
/* generic WinBoard interface, loading its move-generator tables from file */
/***************************************************************************/

     /*****************************************************************/
     /*                      LICENCE NOTIFICATION                     */
     /* Fairy-Max 5.0 is free software, released in the public domain */
     /* so that you have my permission do with it whatever you want,  */
     /* whether it is commercial or not, at your own risk. Those that */
     /* are not comfortable with this, can also use or redistribute   */
     /* it under the GNU Public License or the MIT License.           */
     /* Note, however, that Fairy-Max can easily be configured through*/
     /* its fmax.ini file to play Chess variants that are legally pro-*/
     /* tected by patents, and that to do so would also require per-  */
     /* mission of the holders of such patents. No guarantees are     */
     /* given that Fairy-Max does anything in particular, or that it  */
     /* would not wreck the hardware it runs on, and running it is    */
     /* entirely for your own risk.  H.G,Muller, author of Fairy-Max  */
     /*****************************************************************/

#define MULTIPATH
#define VERSION "5.0b"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <time.h>

#ifndef FAIRYDIR
#define FAIRYDIR "."
#endif
#define INI_FILE FAIRYDIR "/fmax.ini"

#ifdef WIN32 
#    include <windows.h>
#    define CPUtime 1000.*clock
     int Input()
     {  // checks for waiting input in pipe
	static int init; static HANDLE inp; DWORD cnt;
	if(!init) inp = GetStdHandle(STD_INPUT_HANDLE);
	if(!PeekNamedPipe(inp, NULL, 0, NULL, &cnt, NULL)) return 1;
	return cnt;
    }
#else
#    include <sys/time.h>
#    include <sys/times.h>
#    include <sys/ioctl.h>
#    include <unistd.h>
     int GetTickCount() // with thanks to Tord
     {	struct timeval t;
	gettimeofday(&t, NULL);
	return t.tv_sec*1000 + t.tv_usec/1000;
     }
     double CPUtime()
     {  // get CPU time used by process, converted to 'MILLICLOCKS'
	struct tms cpuTimes;
	static int cps = 0;
	if(!cps) cps = sysconf(_SC_CLK_TCK);
	times(&cpuTimes);
	return ((double)(cpuTimes.tms_utime + cpuTimes.tms_stime) * CLOCKS_PER_SEC * 1000)/cps;
     }
     int Input()
     {
	int cnt;
	if(ioctl(0, FIONREAD, &cnt)) return 1;
	return cnt;
     }
#endif

int StartKey;

#define ANALYZE -2
#define EMPTY -1
#define WHITE 0
#define BLACK 16

#define STATE 256

/* The following macros indicate the differences between Fairy-Max and its */
/* dedicated Shatranj derivative ShaMax so that these can now be compiled  */
/* from the same unified source file.                                      */
/* Compile with gcc option -DSHATRANJ to build ShaMax.                     */
#ifdef SHATRANJ
#    define FAC 175
#    define EG  13
#    define NAME "ShaMax"
#    define SHAMAX(x) x
#    define FMAX(x)
#else
#    define FAC 128
#    define EG  10
#    define NAME "Fairy-Max"
#    define SHAMAX(x)
#    define FMAX(x) x
#endif

/* make unique integer from engine move representation */
#define PACK_MOVE 256*K + L + (PromPiece << 16) + (GT<<24);

/* convert intger argument back to engine move representation */
#define UNPACK_MOVE(A) K = (A)>>8 & 255; L = (A) & 255; PromPiece = (A)>>16 & 255; GT = (A)>>24 & 255;

/* Global variables visible to engine. Normally they */
/* would be replaced by the names under which these  */
/* are known to your engine, so that they can be     */
/* manipulated directly by the interface.            */

int Side;
int Move;
int PromPiece;
int Result;
int TimeLeft;
int MovesLeft;
int MaxDepth;
int Post;
int Fifty;
int GameNr;
int Randomize;
int Resign;
char Cambodian[80] = "makruk";
int Threshold = 800;
int drawMoves = 50;
int Score;
int zone, pRank, popup;
int prom, pm, gating, succession, hill;
char piecename[32], piecetype[32], blacktype[32];
char selectedFairy[80];
char *inifile = INI_FILE;
char info[999], hashfile[256];

int Ticks, tlim, Setup, SetupQ;

int GameHistory[1024];
char HistoryBoards[1024][STATE], setupPosition[290];
int GamePtr, HistPtr;
int map[1<<16];

#define RBITS 12
#define W while
#define K(A,B) *(int*)(T+A+S*(B&31))
#define J(A) K(y+A,b[y])-K(x+A,u)-K(H+A,t)
#define SETKEY(Key,A) for(Key=i=0;i<=BE;i++)Key+=K(i+A,b[i])

int U=(1<<23)-1;
struct _ {int K,V;unsigned char X,Y,D,F;} *A;  /* hash table, 16M+8 entries*/

int M=136,S=256,I=8e3,Q,O,K,N,j,R,J,Z,LL,GT,   /* M=0x88                   */
BW,BH,BE,sh,RR,ab,CONS,L,ep,stale,wk,bk,bareK,bareL,score,R2,K2,L2,t3,
pt[513],                                       /* promotion bonus/upgrade  */
lc[32],                                        /* piece location           */
w[16]={0,2,2,-1,7,8,12,23,7,5},                /* relative piece values    */
o[256],
oo[32],                                        /* initial piece setup      */
of[256],
od[16];                                        /* 1st dir. in o[] per piece*/

signed char pl[32],
b[1025],                                       /* board: 16x8+dummy, + PST */
T[8200],                                       /* hash translation table   */
centr[32],
n[]=".*XKNBRQEWFMACHG?x+knbrqewfmachg";        /* piece symbols on printout*/

int pv[10000],*sp=pv; // triangular array
int margin;

int seed=76596595; int Rand() { return (seed = 1103515245*seed + 12345)*150610563>>14; }

void pboard()
{int i;
 i=-1;W(++i<S)printf(" %c",(i&15)==BW&&(i+=15-BW)?10:n[b[i]&31]);
}

TC() { return GetTickCount() - Ticks; }
         
D(k,q,l,e,E,z,n)        /* recursive minimax search, k=moving side, n=depth*/
int k,q,l,e,E,z,n;      /* (q,l)=window, e=current eval. score, E=e.p. sqr.*/
{                       /* e=score, z=prev.dest; J,Z=hashkeys; return score*/
 int j,r,m,v,d,h,i,F,G,P,V,f=J,g=Z,C,s,flag,FF,*ps=sp,kk=S,x,y,X,Y,H,B,rk;
 signed char t,p,u,gt,rg,vf;
 struct _*a=A+(J+(k+S)*E&U);                   /* lookup pos. in hash table*/
 q-=q<e;l-=l<=e;                               /* adj. window: delay bonus */
 d=a->D;m=a->V;X=a->F;Y=a->Y+S-1;              /* resume at stored depth   */
 if(a->K-Z|z&S&&(X=8)||                        /* miss: other pos. or empty*/
  !(m<=q|X&4&&m>=l|X&2))                       /*   or window incompatible */
  d=Y=0;                                       /* start iter. from scratch */
 if(X&1)return 0;                              /* busy-flag set: rep-draw  */
 *sp++=0;                                      /* initialize empty PV      */
 X=a->X;                                       /* start at best-move hint  */
 W(d++<n||d<3||              /*** min depth = 2   iterative deepening loop */
   z&S&&K==I&&(TC()<tlim&d<=MaxDepth||         /* root: deepen upto time   */
   (K=X,L=Y&~S,Score=m,d=3)))                  /* time's up: go do best    */
 {x=B=X;                                       /* start scan at prev. best */
  h=Y&S;                                       /* request try noncastl. 1st*/
  if(a->D<99)a->F=1,a->K=Z;                    /* mark hash entry 'busy'   */
  P=d>2&&l+I?D(16-k,-l,1-l,-e,S,2*S,d-3):I;    /* search null move         */
  m=-P<l|R<5?d-2?-I:e:-P;   /*** prune if > beta  unconsidered:static eval */
  SHAMAX( if(pl[k]<=1&pl[16-k]>1)m=I-1; )      /* bare king loses          */
  ab|=!(N++&4095)&&(tlim>1e8?Input():TC()>t3); /* node count (for timing)  */
  do{u=b[x];                                   /* scan board looking for   */
   if(u&&(u&16)==k)                            /*  own piece (inefficient!)*/
   {r=p=u&15;                                  /* p = piece type (set r>0) */
    if(hill&&w[p]<0&b[769+x])m=I,d=98;         /* King on the hill: we won */
    j=od[p];                                   /* first step vector f.piece*/
    W(r=o[++j])                                /* loop over directions o[] */
    {A:                                        /* resume normal after best */
     flag=h?3:of[j];                           /* move modes (for fairies) */
     y=x;F=FF=G=S;rg=flag>>10&3;vf=32;         /* (x,y)=move, (F,G)=castl.R*/
     if(rg>p)rg=pt[x]&~u&2*u?(vf=0):1;
     do{                                       /* y traverses ray, or:     */
      H=y=h?Y^h:y+r;                           /* sneak in prev. best move */
      if(flag&1<<8)H=y=(y&15)>13?y+BW:(y&15)>=BW?y-BW:y; /* cylinder board */
      if(y<0|y>BE|(y&15)>=BW)break;            /* board edge hit           */
#ifdef MULTIPATH
      if(flag&1<<9)                            /* if multipath move        */
      {t=flag>>RBITS;                          /* get dir. stepped twice   */
       if(b[x+t]){if(b[y-2*t]|b[y-t])break;}else 
       if(b[x+2*t]&&b[y-t])break;              /* test if empty path exists*/
      }
#endif
      m=E<S&&(y<(z&S-1)?E-y<2:y-E<2)&flag?I:m; /* bad castling             */
      if(p<3&flag)H=(y^E)&(E>>9^511)?H:z&S-1;  /* shift capt.sqr. H if e.p.*/
      t=b[H];
      if(flag&1+!t)                            /* mode (capt/nonc) allowed?*/
      {if(t&&(t&16)==k)break;                  /* capture own              */
       i=w[t&15]+((t&192)>>sh);                /* value of capt. piece t   */
       if(i<0){
        if(i<-665)break;                       /* iron piece ends ray scan */
        if(pl[t&31]<2||                        /* K capture, (of last K),  */
        t>>3&kk!=H&kk!=S||(kk=H,i=-i,0))m=I,d=98; /* or duple check: cutoff*/
       }
       if(m>=l&d>1)goto C;                     /* abort on fail high       */
       v=d-1?e:i-p;                            /*** MVV/LVA scoring if d=1**/
       if(d-!t>1)                              /*** all captures if d=2  ***/
       {v=gt=0;G:                              /* retry move with gating   */
        v+=centr[p]*(b[x+513]-b[y+513]);       /* center positional pts.   */
        if(G-S)b[FF]=(rk=b[G])|32,v+=20;       /* castling: put R & score  */
        b[G]=b[H]=0;b[x]=gt;b[y]=u|32;         /* do move, set non-virgin  */
        pl[t&31]-=!!t;                         /* updat victim piece count */
        lc[p+k]=y;                             /* and location             */
        v-=w[p]>0|R<EG?0:20-30*((x-y+1&7)>2);  /*** freeze K in mid-game ***/
        if(p<3)                                /* pawns:                   */
        {v-=9*((b[x-2]!=u)+                    /* structure, undefended    */
               (b[x+2]!=u)                     /*        squares plus bias */
              +(w[b[x^16]&15]<0))              /*** cling to magnetic K ***/
              +(R-76>>2);                      /* end-game Pawn-push bonus */
         b[y]+=V=u&32?pt[y]:0;                 /*upgrade P or convert to Q */
         V>>=sh;                               /* for Shatranj promo to F  */
         i+=V+abs(w[b[y]&15])-w[p];            /* promotion / passer bonus */
        }
        if(z&S){
         if(map[x+S*y]){v=-I;goto S;}          /* skip if excluded move    */
         if(GamePtr<6&Randomize)v+=(Rand()>>10&31)-16;   /* randomize      */
        }
        J+=J(0);Z+=J(4)+G-S;
        SHAMAX( pl[k]-=!!t; )                  /* count pieces per side    */
        v+=e+i;V=m>q?m:q;                      /*** new eval & alpha    ****/
        if(z&S)V=m-margin>q?m-margin:q;        /* multiPV                  */
        C=d-1-(d>5&p>2&!t&!h);                 /* nw depth, reduce non-cpt.*/
        C=R<EG|P-I|d<3||t&&w[p]>0?C:d;         /* extend 1 ply if in-check */
        if(bareK)C=p==bareK&&2&b[769+x]?d+1:d-1; /* corner-leave extension */
        do
         s=C>2|v>V?-D(16-k,-l,-V,-v,/*** futility, recursive eval. of reply */
                                     F,y&255,C):v;
        W(s>q&++C<d); v=s;                     /* no fail:re-srch unreduced*/
        if(v>V&v<l){int *p=sp;
         sp=ps+1;
         W(*sp++=*p++);
         *ps=512*x+y;
        }
        if(z&8*S&&K-I)                         /* move pending: check legal*/
        {if(v+I&&x==K&y==L&gt==GT)             /*   if move found          */
         {Q=-e-i;O=F;LL=L;prom=gt&15;
          if(b[y]-u&15)prom=b[y]-=PromPiece,   /* (under-)promotion:       */
           Q-=abs(w[prom&=15]),Q+=w[prom+PromPiece], /*  correct piece & score & */
                       Z+=PromPiece;           /*  invalidate hash         */
          a->D=99;a->V=0;                      /* lock game in hash as draw*/
          R2-=i/FAC;                           /*** total captd material ***/
          Fifty = t|p<3?0:Fifty+1;
          if(centr[p]>2)bareL=y;               /* remember location bare K */
          sp=ps;
                     return l;}                /*   & not in check, signal */
         v=m;                                  /* (prevent fail-lows on    */
        }                                      /*   K-capt. replies)       */
        SHAMAX( pl[k]+=!!t; ) S:pl[t&31]+=!!t;
        lc[p+k]=x;                             /* restore location of mover*/
        b[G]=rk;b[FF]=b[y]=0;b[x]=u;b[H]=t;    /* undo move,G can be dummy */
       }                                       /*          if non-castling */
       if(z&S&&!ab&K==I&d>2&v>V&v<l){int *p=ps;char X,Y;
        if(Post){
         printf("%2d ",d-2); score=v;
         printf("%6d ", v > I-S ? 100000+I-v : v < S-I ? -100000-I-v : v);
         printf("%8d %10d",(GetTickCount()-Ticks)/10,N);
         while(*p){X=*p>>9;Y=*p++;
         printf(" %c%d%c%d",'a'+(X&15),BH-(X>>4&15)-(BH==10),'a'+(Y&15),BH-(Y>>4&15)-(BH==10));}
         printf("\n");fflush(stdout);
        }GT=gt;K2=x;L2=y;                      /* In root, remember gated  */
       }
       if(v>m)                                 /* new best, update max,best*/
        m=v,X=x,Y=y|S&F;                       /* mark non-double with S   */
       if(gating&&!(u&32)&&p>2&&d-!t>1){       /* virgin non-Pawn: gate    */
        pl[(gt|=k+40)-27]++;                   /* prev. gated back in hand */
        if(m>=l)goto C;                        /* loop skips cutoff :-(    */
        W(++gt<k+43)if(pl[gt-27]){             /* look if more to gate     */
         pl[gt-27]--;v=10;goto G;              /* remove from hand & retry */
       }}
       J=f;Z=g;                                /* restore hash keys        */
       if(ab){a->F&=6;sp=ps;return 0;}         /* unwind search to abort   */
       if(h){h=0;goto A;}                      /* redo after doing old best*/
      }
      s=t&&2&~rg|~t&16^k;v=r^flag>>RBITS;      /* platform & toggled vector*/
      if(flag&15^4|u&vf||                      /* no double or moved before*/
         p>2&!(flag&128)&&                     /* no P & no virgin jump,   */
         ((b[G=x&~15|(r>0)*(BW-1)]^32)<33      /* no virgin R in corner G, */
         ||b[G-r]|b[G-2*r]|b[FF=y+v-r]|b[y+r]) /* no 2 empty sq. next to R */
        )t+=flag&4;                            /* fake capt. for nonsliding*/
      else if(flag&64)t=flag&128?0:t,flag&=63; /* e.p-immune initial step  */
      else F=y+(p<3)*(ep&~u<<8);               /* set e.p. rights          */
      if(s&&flag&8&&!(y=rg&1?y-r:y,t=0)        /* hoppers go to next phase */
         ||!(flag&128)&&!rg--)                 /* zig-zag piece? (w. delay)*/
       r=v,flag^=flag>>4&15;                   /* alternate vector & mode  */
     }W(!t);                                   /* if not capt. continue ray*/
   }}
   if((++x&15)>=BW)x=x>BE?0:x+16&~15;          /* next sqr. of board, wrap */
  }W(x-B);           
C:FMAX( m=m+stale|P==I?m:(X=Y=0); )            /* if stalemate, draw-score */
  if(a->D<99)                                  /* protect game history     */
   a->K=Z,a->V=m,a->D=d,a->X=X,                /* always store in hash tab */
   a->F=4*(m>q)|2*(m<l),a->Y=Y&S?Y+1:0;        /* move, type (bound/exact),*/
 }                                             /*    encoded in X 2,4 bits */
 if(z&4*S)K=X,L=Y&~S;
 sp=ps;
 return m+=m<e;                                /* delayed-loss bonus       */
}


/* Generic main() for Winboard-compatible engine     */
/* (Inspired by TSCP)                                */
/* Author: H.G. Muller                               */

/* The engine is invoked through the following       */
/* subroutines, that can draw on the global vaiables */
/* that are maintained by the interface:             */
/* Side         side to move                         */
/* Move         move input to or output from engine  */
/* PromPiece    requested piece on promotion move    */
/* TimeLeft     ms left to next time control         */
/* MovesLeft    nr of moves to play within TimeLeft  */
/* MaxDepth     search-depth limit in ply            */
/* Post         boolean to invite engine babble      */

/* InitEngine() progran start-up initialization      */
/* InitGame()   initialization to start new game     */
/*              (sets Side, but not time control)    */
/* Think()      think up move from current position  */
/*              (leaves move in Move, can be invalid */
/*               if position is check- or stalemate) */
/* DoMove()     perform the move in Move             */
/*              (togglese Side)                      */
/* ReadMove()   convert input move to engine format  */
/* PrintMove()  print Move on standard output        */
/* Legal()      check Move for legality              */
/* ClearBoard() make board empty                     */
/* PutPiece()   put a piece on the board             */

/* define this to the codes used in your engine,     */
/* if the engine hasn't defined it already.          */

int PrintResult(int s, int mode)
{
        int j, k, cnt=0;

        /* search last 50 states with this stm for third repeat */
        for(j=2; j<=100 && j <= HistPtr; j+=2)
        {
            for(k=0; k<STATE; k++)
                if(HistoryBoards[HistPtr][k] !=
                   HistoryBoards[HistPtr-j&1023][k] )
                   {
                     goto differs;}
            /* is the same, count it */
            if(++cnt > 1) /* third repeat */
            {
                if(mode != EMPTY) printf("1/2-1/2 {Draw by repetition}\n");
                return 1;
            }
          differs: ;
        }
        K=I;ab=0;t3=1e8;
        cnt = D(s,-I,I,Q,O,LL|4*S,3);
#ifdef SHATRANJ
        if(pl[s]==1 && pl[16-s]==1) {
                printf("1/2-1/2 {Insufficient mating material}\n");
                return 4;
        }
        if(pl[s]<=1 && pl[16-s]>1) {
                if (s == BLACK)
                        printf("0-1 {Bare King}\n");
                else
                        printf("1-0 {Bare King}\n");
                return 5;
        }
#else
        if(cnt>-I+1 && K==0 && L==0) {
                printf("1/2-1/2 {Stalemate}\n");
                return 2;
        }
#endif
        if(cnt==-I+1) {
                if (s == WHITE)
                        printf("0-1 {Black mates}\n");
                else {
                        if(succession) { // suppress loss claim if black might be able to replace its King by promotion
                            for(j=0;j<BW;j++)if((b[j+96]&31)==18)return 0;
                        }
                        printf("1-0 {White mates}\n");
                }
                return 3;
        }
        if(Fifty >= 2*drawMoves) {
                if(mode != EMPTY) printf("1/2-1/2 {Draw by fifty move rule}\n");
                return 4;
        }
        return 0;
}


void InitEngine()
{
 N=32*S+7;W(N-->S+3)T[N]=Rand()>>9;
 seed=GetTickCount();
}

void InitGame()
{
 int i,k=0;

 Side = WHITE; Q=0; O=S;
 Fifty = 0; R = 0;
 for(i=0;i<16*BH;i++)b[i]=0;
 for(i=0;i<32;i++)pl[i]=0;
 K=BW;W(K--)
 {b[K]=oo[K+16]+16;b[K+(BH-1)*16]=oo[K];b[K+16*pRank]=18;b[K+(BH-1-pRank)*16]=1; /* initial board setup*/
  pl[oo[K+16]+16]++;pl[oo[K]]++;pl[18]++;pl[1]++;
  if(w[oo[K+16]+16] == -1)pl[oo[K+16]+16]=1;
  if(w[oo[K]] == -1)pl[oo[K]]=1;
  L=BH;W(L--)b[16*L+K+513]=(K-BW/2+hill/2.)*(K-BW/2+hill/2.)+(L-(BH-1)/2.)*(L-(BH-1)/2.),pt[16*L+K]=0; /* center-pts table   */
  pt[K+16]=pt[K+32]=pt[K+(BH-3)*16]=pt[K+(BH-2)*16]=64;pt[K+16*zone]=6-128;pt[K+(BH-1-zone)*16]=5-128; /* promotion bonus & piece upgrade */
  if(pRank == 3) L=oo[K-(w[oo[K]]<0)],pt[K]=L-129,pt[K+(BH-1)*16]=L-130;
 }
 b[769+16*3+BW/2]=b[769+16*4+BW/2]=b[769+16*3+BW/2-1]=b[769+16*4+BW/2-1]=1; /* hill */
 b[769]=b[769+16*(BH-1)]=b[769+16*(BH-1)+BW-1]=b[769+BW-1]=2;            /* corners */
 for(i=0; i<BW; i++) {
  R += abs(w[oo[i]])/FAC + abs(w[oo[i+16]])/FAC;
  Q += abs(w[oo[i]]) - abs(w[oo[i+16]]) + w[1] - w[2];
  if(w[oo[i]] < 0) k = w[oo[i]];
 }
 RR = R2 = R -= 2*(-k/FAC);
 pl[WHITE] = pl[BLACK] = 2*BW; 
 pm = !pl[BLACK+7] && pl[BLACK+9] && pl[WHITE+7] ? 2 : 0; // Unlike white, black has no 'Q', so promote to 9, which he does have.
 K=BW;W(K--)pt[K+(BH-1)*16]+=pm; // alter black promotion choice
 if(gating) pl[14] = pl[15] = pl[30] = pl[31] = 1, R2 = R += 2*(w[9]/FAC + w[10]/FAC);
}

void CopyBoard(int s)
{
        int i, j;

        /* copy game representation of engine to HistoryBoard */
        /* don't forget castling rights and e.p. state!       */
        for(i=0; i<BH; i++)
        for(j=0; j<BW; j++)                 /* board squares  */
            HistoryBoards[s][BW*i+j] = b[16*i+j]|64*(16*i+j==O);
}
                                         
void PrintVariants(int combo)
{
        int count=0, total=0, c=EOF+1; char buf[80];
        FILE *f;

        f = fopen(inifile, "r");
        if(f==NULL) return;

        /* search for game names in definition file */
        do {
           while(fscanf(f, "Game: %s", buf) != 1 && c != EOF) 
               while((c = fgetc(f)) != EOF && c != '\n');
           if(c == EOF) break;
           total++;
           if(*buf < '1' || *buf >'9' && *buf < 'a' || combo == (strstr(buf, "fairy/") != buf)) continue;
           if(combo && count == 0) strcpy(selectedFairy, buf);
           if(count++) printf(combo ? " /// " : ",");
           printf("%s", combo ? buf+6 : buf);
        } while(c != EOF);

        fclose(f);
        if(!combo && total != count) printf("%sfairy", count ? "," : "");
}

void PrintOptions()
{
	printf("feature option=\"Resign -check %d\"\n", Resign);
	printf("feature option=\"Resign Threshold -spin %d 200 1200\"\n", Threshold);
	printf("feature option=\"Claim draw after -spin %d 0 200\"\n", drawMoves);
	printf("feature option=\"Ini File -file %s\"\n", inifile);
	printf("feature option=\"Multi-PV Margin -spin %d 0 1000\"\n", margin);
	printf("feature option=\"Variant fairy selects -combo "); PrintVariants(1); printf("\"\n");
	printf("feature option=\"Makruk rules -combo makruk /// Cambodian /// Ai-wok\"\n");
	printf("feature option=\"Dummy Slider Example -slider 20 0 100\"\n");
	printf("feature option=\"Dummy String Example -string happy birthday!\"\n");
	printf("feature option=\"Dummy Path Example -path .\"\n");
	printf("feature option=\"Automatic persistent-hash dialog -check %d\"\n", popup);
	printf("feature option=\"Info -button\"\n");
	printf("feature option=\"Save in hash file -button\"\n");
	printf("feature option=\"Clear Hash -button\"\n");
	printf("feature done=1\n");
}
                                         
void LoadHash(char *dir, char *name)
{       // read persistent-hash file into hash table as protected entries that always cause cutoff
        FILE *f;
        snprintf(hashfile, 256, "%s/%s.hash", dir, name);
        if(f = fopen(hashfile, "r")) {
            while(fscanf(f, "%x:%x=%d", &J, &Z, &j) == 3) J&=U, A[J].V=j, A[J].K=Z, A[J].D=99, A[J].F=6;
            fclose(f);
        }
}

void LoadGame(char *name)
{
        int i, j, ptc=0, count=0, step2; char c, buf[80], pieceToChar[200], parent[80];
        static int currentVariant;
        FILE *f;

        f = fopen(inifile, "r");
        if(f==NULL)
        {   printf("telluser piece-description file '%s'  not found\n", inifile);
            exit(0);
        }
        if(fscanf(f, "version 4.8(%c)", &c)!=1 || c != 'w')
        { printf("telluser incompatible fmax.ini file\n"); exit(0); }

        gating = succession = 0;
        if(name != NULL)
        {  /* search for game name in definition file */
           if(!strcmp(name, "makruk")) name = Cambodian; else
           if(!strcmp(name, "fairy")) name = selectedFairy;
           gating = !strcmp(name, "seirawan");
           while((ptc=fscanf(f, "Game: %s # %s # %s", buf, pieceToChar, parent))==0 || strcmp(name, buf) ) {
               char *p = info;
               while((c = fgetc(f)) != EOF && c != '\n') *p++ = c;
               if(*info == '/') *p = 0; else *info = 0; // remember last line before Game if it was comment
               count++;
               if(c == EOF) {
                   printf("telluser variant %s not supported\n", name);
                   fclose(f);
                   return; /* keep old settings */
               }
           }
           currentVariant = count;
        }

        /* We have found variant, or if none specified, are at beginning of file */
        if(fscanf(f, "%dx%d", &BW, &BH)!=2 || BW>14 || BH>16)
        { printf("telluser unsupported board size %dx%d\n",BW,BH); exit(0); }
        BE = (BH-1)*16 + BW-1; CONS = 799 + 16*(BH-8);   // highest valid square number and move-conversion constant
        i = 1; fscanf(f, "=%d", &i); zone = i - 1;       // new method to indicate deviant zone depth

        for(i=0; i<BW; i++) fscanf(f, "%d", oo+i);
        for(i=0; i<BW; i++) fscanf(f, "%d", oo+i+16);
        for(i= 0; i<=U; i++)
            A[i].K = A[i].D = A[i].X = A[i].Y = A[i].F = 0; /* clear hash */
        for(i=0; i<32; i++) piecetype[i] = blacktype[i] = 0;

        i=0; j=-1; c=0; ep=1<<20; stale=I; bk=1; bareK=0; bareL=-1; step2 = 666;
        while(fscanf(f, "%d,%x,%d", o+j, of+j, &step2)>=2 ||
                                      fscanf(f,"%c:%d",&c, w+i+1)==2)
        {   if(c)
            { od[++i]=j; centr[i] = c>='a';
              blacktype[c&31]=i; piecename[i]=c&31;
              if(piecetype[c&31]==0) piecetype[c&31]=i; // only first
              succession |= w[i] < -4;         // expendable royalty; assume we can promote to it
              if(w[i]<0) wk=bk, bk=i;  // remember royals
            }
            if(step2 != 666) of[j] += (step2 ^ o[j]) << RBITS, step2 = 666; // compute toggle-vector from 3rd move parameter
            j++; o[j]=0;
            /* printf("# c='%c' i=%d od[i]=%d j=%d (%3d,%8x)\n",c?c:' ',i,od[i],j,o[j-1],of[j-1]); /**/
            c=0; if(i>15 || j>255) break;
        }

	if(BH == 10 && o[0] == -16 && of[0] & 0xC00) ep += 16<<9; // pawn with triple-push
	sh = w[7] < 250 ? 3 : 0; hill = (w[3] == -2); stale -= (w[9] == -2); pRank = (zone ? zone : 1);
        if(zone < 0) pRank = -1-zone, zone = 0; // negative =N suffix is kludge for configuring Pawn rank
        if(ptc > 1) { // setup board in GUI, by sending it pieceToCharTable and FEN
            if(ptc == 2) printf("setup (%s) ", pieceToChar);
            else printf("setup (%s) %dx%d+0_%s ", pieceToChar, BW, BH, parent);
            for(i=0; i<BW; i++) printf("%c", piecename[oo[i+16]]+'`'); printf("/");
            for(i=1; i<pRank; i++) printf("8/");
            for(i=0; i<BW; i++) printf("%c", piecename[2]+'`'); printf("/");
            for(i=1+pRank; i<BH-1-pRank; i++) printf("%d/", BW);
            for(i=0; i<BW; i++) printf("%c", piecename[1]+'@'); printf("/");
            for(i=1; i<pRank; i++) printf("8/");
            for(i=0; i<BW; i++) printf("%c", oo[i] ? piecename[oo[i]]+'@' : '1'); printf(" w KQkq - 0 1\n");
        }
	while(fscanf(f, " # %[^\n]", pieceToChar)) printf("piece %s\n", pieceToChar);
        fclose(f);
        LoadHash(FAIRYDIR, name); LoadHash(".", name); // initialize persistent hash
}

int main(int argc, char **argv)
{
        int Computer, MaxTime, MaxMoves, TimeInc, sec, i;
        char line[256], command[256], c, ff, ft;
        int m, nr, rf, rt;
        double cpuT;

        if(argc > 1 && !strcmp(argv[1], "-v")) argc++, argv--, printf("%s\n", VERSION), exit(0);

        if(argc>1 && sscanf(argv[1], "%d", &m)==1)
        { U = (1<<m)-1; argc--; argv++; }
        A = (struct _ *) calloc(U+1, sizeof(struct _));
        if(argc>1) inifile = argv[1];

	signal(SIGINT, SIG_IGN);
	setvbuf(stdin, NULL, _IONBF, 0); // suppress input buffering
        printf("tellics say     " NAME " " VERSION "\n");
        printf("tellics say     by H.G. Muller\n");
        InitEngine();
        LoadGame(NULL);
        InitGame();
        Computer = EMPTY;
        MaxTime  = 10000;  /* 10 sec */
        MaxDepth = 30;     /* maximum depth of your search */

        for (;;) {
		fflush(stdout);
                PromPiece = 0; /* Always promote to Queen ourselves */
                for(N=K=0;K<S;K++)N+=b[K]?b[K]&16?S:1:0; /* count pieces for detecting bare King */
                if(w[wk]<0&w[bk]<0){if(N<2*S)bareK=bk;if(!(N&S-2))bareK=wk;}
                R = R2 - 2*abs(Q)/(3*FAC); if(R < 0) R=0;  /* treat strongly unbalanced as if later game phase */
                if(bareK)centr[bareK]=1+Fifty/10,R=4;
                SETKEY(J,0); SETKEY(Z,4); /* absolutize key, so it can be used for persistent hash */
                K=(bareL&15);L=bareL>>4;
                if(!K|K==BW-1&&!L|L==BH-1&&b[513+bareL])
                 for(i=0;i<BH;i++)for(m=0;m<BW;m++)b[513+16*i+m]=abs(abs(i-L)-abs(m-K));
                if(hill) centr[3] = R>20 ? 1 : 22-R;
                Ticks = GetTickCount();
                if (Side == Computer) {
                        /* think up & do move, measure time used  */
                        /* it is the responsibility of the engine */
                        /* to control its search time based on    */
                        /* MovesLeft, TimeLeft, MaxMoves, TimeInc */
                        /* Next 'MovesLeft' moves have to be done */
                        /* within TimeLeft+(MovesLeft-1)*TimeInc  */
                        /* If MovesLeft<0 all remaining moves of  */
                        /* the game have to be done in this time. */
                        /* If MaxMoves=1 any leftover time is lost*/
                        cpuT = CPUtime(); printf("# times @ %u\n", Ticks);
                        m = MovesLeft<=0 ? 40 : MovesLeft;
                        t3 = (TimeLeft+(m-1)*TimeInc)/(m+2);
                        tlim = (0.6-0.06*(BW-8))*t3; t3 *= 3; t3 -= 32;
                        if(tlim>TimeLeft/15) tlim = TimeLeft/15;
printf("# %d+%d pieces, centr = (%d,%d) R=%d\n", N&63, N>>8, centr[wk], centr[bk], R);
                        if(bareK|RR>4^R>4) // with bare King or after switching on or off null move
                            for(i=0;i<=U;i++)if(A[i].D<99&&abs(A[i].V)<I-S)A[i].K=0; // clear hash
                        N=ab=0;K=I;RR=R;
                        if (D(Side,-I,I,Q,O,LL|9*S,3)==I || ab && 
                            (K=K2,L=L2,ab=0,t3=1e8,D(Side,-I,I,Q,O,LL|9*S,3))==I) { // make move for aborted search
                            Side ^= BLACK^WHITE;
                            m = GetTickCount() - Ticks;
                            printf("# times @ %u: real=%d cpu=%1.0f\n", m + Ticks, m,
                                      (CPUtime() - cpuT)/CLOCKS_PER_SEC);
printf("# promo = %d (%c) GT = %d\n", prom, piecename[prom]+'`', GT); 
                            printf("move ");
                            printf("%c%d%c%d",'a'+(K&15),BH-(K>>4)-(BH==10),
                                              'a'+(L&15),BH-(L>>4)-(BH==10));
			    if(prom)printf("%c",piecename[prom]+'a'-1);
                            printf("\n");

                            /* time-control accounting */
                            TimeLeft -= m;
                            TimeLeft += TimeInc;
                            if(--MovesLeft == 0) {
                                MovesLeft = MaxMoves;
                                if(MaxMoves == 1)
                                     TimeLeft  = MaxTime;
                                else TimeLeft += MaxTime;
                            }

                            GameHistory[GamePtr++] = PACK_MOVE;
                            CopyBoard(HistPtr=HistPtr+1&1023);
			    if(Resign && Score <= -Threshold) { 
				printf("resign\n"); Computer=EMPTY;
                            } else if(PrintResult(Side, Computer))
                                Computer = EMPTY;
                        } else {
                            if(!PrintResult(Side, Computer))
                                printf("resign { refuses own move }\n");
                            Computer = EMPTY;
                        }
                        continue;
		}
		if(Computer == ANALYZE) {
			if(popup-- == 1) popup++, printf("askuser remember Save score in hash file (OK/Cancel)?\n");
                        N=ab=0;K=I;tlim=t3=1e9;
			D(Side,-I,I,Q,O,LL|S,3);
		}
		if (!fgets(line, 256, stdin))
			return 1;
		if (line[0] == '\n')
			continue;
		sscanf(line, "%s", command);
		if (!strcmp(command, "xboard"))
			continue;
                if (!strcmp(command, "protover")) {
                        printf("feature myname=\"" NAME " " VERSION "\"\n");
                        printf("feature memory=1 exclude=1\n");
                        printf("feature setboard=0 xedit=1 ping=1 done=0\n");
                        printf("feature variants=\"");
                        PrintVariants(0);
                        printf("\"\n");
			PrintOptions();
                        continue;
                }
                if (!strcmp(command, "ping")) { int nr=0;
                        sscanf(line, "ping %d", &nr);
                        printf("pong %d\n", nr);
			continue;
                }
                if (!strcmp(command, "p")) {
                        pboard();
			continue;
                }
                if (!strcmp(command, "memory")) {
                        int mem, mask;
			sscanf(line+6, "%d", &mem); mem = (mem*1024*1024)/12; // max nr of hash entries
			mask = 0x3FFFFFF; while(mask > mem) mask >>= 1;
			if(mask != U) {
			    free(A); U = mask;
			    A = (struct _ *) calloc(U+1, sizeof(struct _));
			}
			continue;
                }
#		define CLEAR(X) for(i=0; i<S*S; i++) map[i] = X
		if (!strcmp(command+2, "clude")) { // include / exclude
			char *c=line+8, K=c[0]-16*c[1]+CONS, L=c[2]-16*c[3]+CONS, r = *command - 'i';
			if(!strcmp(line+8, "all\n")) CLEAR(r);
			else map[K+S*L] = r;
			continue;
		}
		CLEAR(0);
		if (!strcmp(command, "new")) {
                        /* start new game */
                        LoadGame("normal");
                        InitGame();
                        GamePtr   = Setup = 0;
                        GameNr++;
                        HistPtr   = 0;
                        Computer  = BLACK;
                        TimeLeft  = MaxTime;
                        MovesLeft = MaxMoves;
                        Randomize = 0;
                        for(nr=0; nr<1024; nr++)
                            for(m=0; m<STATE; m++)
                                HistoryBoards[nr][m] = 0;
			continue;
		}
		if (!strcmp(command, "quit"))
                        /* exit engine */
			return 0;
		if (!strcmp(command, "analyze")) {
                        /* computer plays neither */
                        Computer = ANALYZE; Randomize *= 2;
			continue;
		}
		if (!strcmp(command, "exit") ||
		    !strcmp(command, "force")) {
                        /* computer plays neither */
                        Computer = EMPTY; Randomize = (Randomize > 0);
			continue;
		}
		if (!strcmp(command, "white")) {
                        /* set white to move in current position */
                        if(Side == BLACK) Q = -Q;
                        Side     = WHITE;
                        Computer = BLACK;
			continue;
		}
		if (!strcmp(command, "black")) {
                        /* set blck to move in current position */
                        if(Side == WHITE) Q = -Q;
                        Side     = BLACK;
                        Computer = WHITE;
			continue;
		}
		if (!strcmp(command, "st")) {
                        /* move-on-the-bell mode     */
                        /* indicated by MaxMoves = 1 */
                        sscanf(line, "st %d", &MaxTime);
                        MovesLeft = MaxMoves = 1;
                        TimeLeft  = MaxTime *= 1000;
                        TimeInc   = 0;
			continue;
		}
		if (!strcmp(command, "sd")) {
                        /* set depth limit (remains in force */
                        /* until next 'sd n' command)        */
                        sscanf(line, "sd %d", &MaxDepth);
                        MaxDepth += 2; /* QS depth */
			continue;
		}
                if (!strcmp(command, "level")) {
                        /* normal or blitz time control */
                        sec = 0;
                        if(sscanf(line, "level %d %d %d",
                                 &MaxMoves, &MaxTime, &TimeInc)!=3 &&
                           sscanf(line, "level %d %d:%d %d",
                                 &MaxMoves, &MaxTime, &sec, &TimeInc)!=4)
                             continue;
                        MovesLeft = MaxMoves;
                        TimeLeft  = MaxTime = 60000*MaxTime + 1000*sec;
                        TimeInc  *= 1000;
                        continue;
                }
		if (!strcmp(command, "time")) {
                        /* set time left on clock */
                        sscanf(line, "time %d", &TimeLeft);
                        TimeLeft  *= 10; /* centi-sec to ms */
			continue;
		}
		if (!strcmp(command, "otim")) {
                        /* opponent's time (not kept, so ignore) */
			continue;
		}
		if (!strcmp(command, "easy")) {
			continue;
		}
		if (!strcmp(command, "hard")) {
			continue;
		}
		if (!strcmp(command, "accepted")) {
			continue;
		}
		if (!strcmp(command, "rejected")) {
			continue;
		}
		if (!strcmp(command, "random")) {
			Randomize = !Randomize;
			continue;
		}
		if (!strcmp(command, "remember")) {
			FILE *f = fopen(hashfile, "a"); // add current position to persistent hash
			sscanf(line+8, "%d", &score);   // user can overrule score
			if(f) fprintf(f, "%08x:%08x=%d\n", J+(S+Side)*O, Z, score), fclose(f);
			popup = 2; // suppresses repeat of popup on restart of analysis search
			continue;
               }
		if (!strcmp(command, "option")) {
			int i; static char filename[80];
			if(sscanf(line+7, "Resign=%d", &Resign) == 1) continue;
			if(sscanf(line+7, "Resign Threshold=%d", &Threshold) == 1) continue;
			if(sscanf(line+7, "Ini File=%s", filename) == 1) {
				inifile = filename; continue;
			}
			if(sscanf(line+7, "Clear Hash%c", &c) == 1) for(i=0; i<=U; i++) A->K = 0;
			if(sscanf(line+7, "Info%c", &c) == 1) printf("telluser %s\n", info+3);
			if(sscanf(line+7, "MultiVariation Margin=%d", &margin) == 1) continue;
			if(sscanf(line+7, "Variant fairy selects=%s", selectedFairy+6) == 1) continue;
			if(sscanf(line+7, "Makruk rules=%s", Cambodian) == 1) continue;
			if(sscanf(line+7, "Claim draw after=%d", &drawMoves) == 1) continue;
			if(sscanf(line+7, "Automatic persistent-hash dialog=%d", &popup) == 1) continue;
			if(sscanf(line+7, "Save in hash file%c", &c) == 1 && Computer == ANALYZE) {
			     FILE *f = fopen(hashfile, "a"); // add current position to persistent hash
			     if(f) fprintf(f, "%08x:%08x=%d\n", J+(S+Side)*O, Z, score), fclose(f);
			}
			continue;
		}
		if (!strcmp(command, "go")) {
                        /* set computer to play current side to move */
                        Computer = Side;
                        MovesLeft = -(GamePtr+(Side==WHITE)>>1);
                        while(MaxMoves>0 && MovesLeft<=0)
                            MovesLeft += MaxMoves;
			continue;
		}
		if (!strcmp(command, "hint")) {
                        Ticks = GetTickCount(); tlim = 1000; ab = 0; t3 = 1e8;
                        D(Side,-I,I,Q,O,LL|4*S,6);
                        if (K==0 && L==0)
				continue;
                        printf("Hint: ");
                        printf("%c%d%c%d",'a'+(K&15),BH-(K>>4)-(BH==10),
                                          'a'+(L&15),BH-(L>>4)-(BH==10));
                        printf("\n");
			continue;
		}
                if (!strcmp(command, "undo")   && (nr=1) ||
                    !strcmp(command, "remove") && (nr=2)   ) {
                        /* 'take back' moves by replaying game */
                        /* from history until desired ply      */
                        if (GamePtr - nr < 0)
				continue;
                        GamePtr -= nr;
                        HistPtr -= nr;   /* erase history boards */
                        while(nr-- > 0)  
                            for(m=0; m<STATE; m++)
                                HistoryBoards[HistPtr+nr+1&1023][m] = 0;
                        InitGame();
			if(Setup) {
			    for(i=0; i<S; i++) b[i] = setupPosition[i];
			    for(i=0; i<32; i++) pl[i] = setupPosition[i+S+2];
			    Side = setupPosition[S]; Q = SetupQ;
			    R2 = setupPosition[S+1];
			}
			for(i=0; i<=U; i++) if(A[i].D == 99) A[i].D = A[i].K = 0; // clear game history from hash table
                        for(nr=0; nr<GamePtr; nr++) {
                            UNPACK_MOVE(GameHistory[nr]); SETKEY(J,0); SETKEY(Z,4);
                            ab=0;t3=1e8;D(Side,-I,I,Q,O,LL|9*S,3);
                            Side ^= BLACK^WHITE;
                        }
			continue;
		}
		if (!strcmp(command, "post")) {
                        Post = 1;
			continue;
		}
		if (!strcmp(command, "nopost")) {
                        Post = 0;
			continue;
		}
		if (!strcmp(command, "variant")) {
                        sscanf(line, "variant %s", command);
                        LoadGame(command);
                        InitGame(); Setup = 0;
			continue;
		}
                if (!strcmp(command, "edit")) {
                        int color = WHITE, p, r;

                        while(fgets(line, 256, stdin)) {
                                m = line[0];
                                if(m=='.') break;
                                if(m=='#') {
                                        for(i=0; i<S; i++) b[i]=0;
                                        for(i=0; i<32; i++) pl[i]=0;
                                        Q=0; R=0; O=S;
                                        pl[WHITE]=pl[BLACK]=0;
                                        continue;
                                }
                                if(m=='c') {
                                        color = WHITE+BLACK - color;
                                        Q = -Q;
                                        continue;
                                }
                                if( m >= 'A' && m <= 'Z' && piecetype[m&31]) {
                                    p = (color == WHITE ? piecetype : blacktype)[line[0]&31];
                                    if(line[1] == '@') { // stuff holdings
                                        pl[color+p+5] = m = line[2] - '0';
                                        pl[BLACK+WHITE-color]+=m;pl[p+color]+=m;
                                        Q+=m*w[p]; R+=m*(w[p]/FAC);
                                        continue;
                                    } else
                                    if(line[1] >= 'a' && line[1] <= 'a'+BW-1
                                    && line[2] >= '0' && line[2] <= '0'+BH) {
                                        line[2] = '0' + atoi(line + 2) + (BH==10); // allow 2-digit rank
                                        m = line[1]-16*line[2]+CONS; r = m & 0xF0;
                                        switch(p)
                                        {
                                        case 1:
                                        case 2:
                                            if(color==WHITE)
                                                 b[m]=r==0x10?161:r==0x20?97:r>=16*(BH-1-pRank)?1:33,
                                                 Q+=w[1]+(r==0x10?128:r==0x20?64:0);
                                            else b[m]=r==16*(BH-2)?178:r==16*(BH-3)?114:r<=0x10*pRank?18:50,
                                                 Q+=w[2]+(r==16*(BH-2)?128:r==16*(BH-3)?64:0);
                                            break;
                                        default:
                                            b[m]=p+color+32; // assume non-virgin
					    if(color==BLACK && m<0x10 && p==oo[m+16] || // but make virgin on original square
                                               color==WHITE && m>=16*(BH-1) && p==oo[m-16*(BH-1)]) b[m] -= 32;
                                            if(w[p]<0) { // Royal piece on original square: virgin
                                                Q-=w[p]; // assume value was flipped to indicate royalty
                                                if(pl[p+color])R-=w[p]/FAC; // capturable King, add to material
					    } else { Q+=w[p]; R+=w[p]/FAC; }
					case 0: // undefined piece, ignore
                                            break;
                                        }
                                        pl[BLACK+WHITE-color]++;pl[p+color]++;
                                        if(w[p+color] == -1)pl[p+color]=1; // fake we have one if value = -1, to thwart extinction condition
                                        continue;
                                    }
                                }
                        }
                        if(Side != color) Q = -Q;
			GamePtr = HistPtr = 0; Setup = 1; SetupQ = Q; // start anew
			for(i=0; i<S; i++) setupPosition[i] = b[i]; // remember position
			setupPosition[S] = Side;
			setupPosition[S+1] = RR = R2 = R;
			for(i=0; i<32; i++) setupPosition[i+S+2] = pl[i];
			Computer = EMPTY; // after edit: force mode!
			continue;
		}
                /* command not recognized, assume input move */
                GT = 0; m = (sscanf(line, "%c%d%c%d%c", &ff, &rf, &ft, &rt, &c) < 5);
                if(BH==10)rf++,rt++;
                if(c != '\n') GT = (Side == WHITE ? piecetype : blacktype)[c&31];
                K=ff-16*(rf+'0')+CONS;L=ft-16*(rt+'0')+CONS;
                if(GT) PromPiece = (pt[L]&15) + 1 + (Side == BLACK) - GT, GT |= 32 + Side;
                if(w[GT&15] == -1 || w[GT&15]%10 == 3) L = S; // spoil move for promotion to King (or when marked non-promoting)
                if(pRank == 3 && PromPiece) L = S;            // no promotion choice, spoil move if not default piece
                if (m & line[1] != '@')
                        /* doesn't have move syntax */
			printf("Error (unknown command): %s\n", command);
                else { int i=-1;
                    if(b[L] && (b[L]&16) == Side && w[b[L]&15] < 0) // capture own King: castling
                    { i=K; K = L; L = i>L ? i-1 : i+2; }
		    if(w[GT&15] < -1) pl[GT&31]++, J+=89729; // promotion to royal piece
                    if((b[K]&15) < 3) GT = 0; // Pawn => true promotion rather than gating
                    ab=0;t3=1e8;
                    if(D(Side,-I,I,Q,O,LL|9*S,3)!=I) {
                        /* did have move syntax, but illegal move */
                        printf("Illegal move:%s\n", line);
                    } else {  /* legal move, perform it */
                        if(i >= 0) b[i]=b[K],b[K]=0; // reverse Seirawan gating
                        GameHistory[GamePtr++] = PACK_MOVE;
                        Side ^= BLACK^WHITE;
                        CopyBoard(HistPtr=HistPtr+1&1023);
                        if(PrintResult(Side, Computer) && Computer != ANALYZE) Computer = EMPTY;
		    }

		}
	}
	return 0;
}
