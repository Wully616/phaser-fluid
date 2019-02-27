//https://mikeash.com/pyblog/fluid-simulation-for-dummies.html
//http://www.cs.northwestern.edu/~sco590/cs140/assignment3.html
//https://github.com/dionyziz/wave-experiment/blob/master/main.coffee
//http://curran.github.io/HTML5Examples/
//https://www.thanassis.space/wavePhysics.html
//http://cowboyprogramming.com/2008/04/01/practical-fluid-mechanics/
import Vector from '../../util/Vector2D.js'
export default class FluidGrid extends Phaser.GameObjects.GameObject {

    constructor(scene, x, y, key, gridWidth,gridHeight,blockSize, diffusion, viscosity, dt) {
        super(scene, x, y, key);
        this.scene = scene;
        this.scene.add.existing(this);

        this.gridWidth = gridWidth;
        this.gridHeight = gridHeight;
        this.blockSize = blockSize;
        //a debug grid
        //this.createDebugGrid(this.scene,this.gridWidth,this.gridHeight,this.blockSize);
        
        this.iter = 20;
        //size of the grid
        this.size = gridWidth/blockSize;
        this.gridsize = this.size*this.size;
        this.N = this.size;
        //the timestamp, diffusion and viscosity of this fluidgrid
        this.dt = dt;
        this.diff = diffusion;
        this.visc = viscosity;
        this.sR = new Array(this.gridsize).fill(0);
        this.sG = new Array(this.gridsize).fill(0);
        this.sB = new Array(this.gridsize).fill(0);
        this.densityR = new Array(this.gridsize).fill(0);
        this.densityG = new Array(this.gridsize).fill(0);
        this.densityB = new Array(this.gridsize).fill(0);
        this.Vx = new Array(this.gridsize).fill(0);
        this.Vx0 = new Array(this.gridsize).fill(0);
        this.Vy = new Array(this.gridsize).fill(0);
        this.Vy0 = new Array(this.gridsize).fill(0);

        this.solids = new Array(this.gridsize).fill(false);
        this.tmpar = []; //three temp arrays
        for(var ar = 0 ; ar < 3; ar ++){
            this.tmpar[ar] = new Array(this.gridsize).fill(0);
        }
        this.graphics = this.scene.add.graphics();

        this.wrap = true; // wrap grid bounds
        this.ixList = {};

        for (var i = -10; i < this.N +10 ; i++) {
            this.ixList[i] = {};
            for (var j = -10; j < this.N+10; j++) { 
                
                this.ixList[i][j] = {
                    up: this.IXWrap(i,j - 1),
                    down: this.IXWrap(i,j + 1),
                    left: this.IXWrap(i-1 ,j),
                    right: this.IXWrap(i+1,j),
                    ix: this.IXWrap(i,j)
                }
            }
        }

        //console.log(this.ixList);

    }

    createDebugGrid(scene,width,height,blockSize){
        scene.add.grid(width/2,height/2,width, height,blockSize,blockSize,0,255,200);
        
    }


    fluidStep(){
        var N        = this.size;
        var visc     = this.visc;
        var diff     = this.diff;
        var dt       = this.dt;
        var Vx      = this.Vx;
        var Vy      = this.Vy;
        var Vx0     = this.Vx0;
        var Vy0     = this.Vy0;
        var sR       = this.sR;
        var sG       = this.sG;
        var sB       = this.sB;
        var densityR = this.densityR;
        var densityG = this.densityG;
        var densityB = this.densityB;
            
        //this.diffuseVec(1, Vx0, Vx, 2, Vy0, Vy, visc, dt, this.iter, N );
        this.diffuseArr([Vx0,Vy0 ], [Vx, Vy], visc, dt, this.iter, N );
        this.project(Vx0, Vy0, Vx, Vy, this.iter, N);
        this.edgeVelocities(Vx,Vy, N);
        this.advectArr([Vx,Vy],[Vx0,Vy0], Vx0, Vy0,  dt, N); 
        this.project(Vx, Vy, Vx0, Vy0, this.iter, N);
        this.edgeVelocities(Vx,Vy, N);
        //process R G B channels of density
        this.diffuseArr([sR,sG,sB], [densityR,densityG,densityB], diff, dt, this.iter, N);
        this.advectArr([densityR,densityG,densityB], [sR,sG,sB], Vx, Vy, dt, N);

        
  
    }


    diffuseArr(x, x0, diff, dt, iter, N){
        var a = dt * diff * (N - 2) * (N - 2);
        this.lin_solveArr(x, x0, a, 1 + 4 * a, iter, N);
    }

    lin_solve(x, x0,  a,  c,  iter,  N){
        var cRecip = 1.0 / c;

        for (var k = 0; k < iter; k++) {     
            for (var j = 0; j < N ; j++) {
                for (var i = 0; i < N; i++) {
                    var idx = this.IXObj(i,j);                    
                    this.tmpar[0][ idx.ix ] = this.calcNeighborVal(idx,x,x0,a, cRecip, N);                    
                }
            }
            for(var i = 0; i < N*N; i++){
                x[i] = this.tmpar[0][i];
            }

        }
       

    }
    lin_solveArr(x, x0,  a,  c,  iter,  N){
        var cRecip = 1.0 / c;
        for (var k = 0; k < iter; k++) {     
            for (var j = 0; j < N ; j++) {
                for (var i = 0; i < N; i++) {
                    var idx = this.IXObj(i,j);
                    for(var ar = 0 ; ar < x.length; ar ++){
                        this.tmpar[ar][ idx.ix ] = this.calcNeighborVal(idx ,x[ar],x0[ar],a, cRecip, N);
                    }
                }
            }      
            for(var i = 0; i < N*N; i++){
                for(var ar = 0 ; ar < x.length; ar ++){
                    x[ar][i] = this.tmpar[ar][i];
                }          
            }
            
        }
    }

    calcNeighborVal(i ,x,x0,a, cRecip, N){
        var ix = i.ix;
        var right = i.right;
        var left = i.left;
        var up = i.up;
        var down = i.down;
        if(!this.solids[ix]){

            //create bitshift to see if walls are up down left right,
            //then do if statements to calculate appropriate diffusion based on           
            return (
                    x0[ix] + a *
                    (   this.solids[right]   ? x[ix] : x[ right ]  + //+ 1 on x axis
                        this.solids[left]    ? x[ix] : x[ left ]  + //-1 on x axis
                        this.solids[down]    ? x[ix] : x[ down ]  + /// +1 on y axis
                        this.solids[up]      ? x[ix] : x[ up ] // -1 on y axis
                    )
                ) * cRecip;
        } else {
            return 0; //blocks should not contain densities
        }
            
    }

    project(velocX, velocY, p, div, iter, N)
    {
        var nRecip = 1 / N;   
        for (var j = 0; j < N ; j++) {
            for (var i = 0; i < N; i++) {
                var idx = this.IXObj(i, j);
                div[ idx.ix ] = -0.5*(
                            (  velocX[ idx.right] )
                            -( velocX[ idx.left  ] )
                            +( velocY[ idx.down  ])
                            -( velocY[ idx.up  ] )
                ) * nRecip;
                p[ idx.ix ] = 0;
            }
        }
    
        this.lin_solve( p, div, 1, 6, iter, N);
           
        for (var j = 0; j < N ; j++) {
            for (var i = 0; i < N; i++) {
                var idx = this.IXObj(i, j);
                velocX[ idx.ix ] -= 0.5 * (  p[ idx.right ]- p[ idx.left ]) * N;
                velocY[ idx.ix ] -= 0.5 * (  p[ idx.down ]- p[ idx.up ]) * N;
            }
        }

    }

    advectArr(d, d0, u, v, dt, N)
    {
        var i0, i1, j0, j1;
               
        var s0, s1, t0, t1;
        var  x, y;
        
        var i, j;
        var dt0 = dt*N;
        var tleft, t, tnext,  vx, vy;

        for(j = 0; j < N; j++) { 
            for(i = 0; i < N; i++) {
                var ix = this.IX(i, j);
                //if(this.solids[ix]){continue;}

                tleft=dt0;
                x=i;y=j;		
        
                while(tleft>0.0000001) {
        
                    //enforce boundry contraints
                    if (x<0.5) x=0.5; if (x>N+0.5) x=N+0.5; 
                    if (y<0.5) y=0.5; if (y>N+0.5) y=N+0.5; 
                    
                    var c = this.checkBounds(x,y);
                    i0=c[0]; i1=i0+1;
                    j0=c[1]; j1=j0+1;
                    s1 = x-i0; s0 = 1-s1; t1 = y-j0; t0 = 1-t1;
        
                    vx = -(s0*(t0*u[this.IX(i0,j0)]+t1*u[this.IX(i0,j1)])+
                           s1*(t0*u[this.IX(i1,j0)]+t1*u[this.IX(i1,j1)]));
                            
                    vy = -(s0*(t0*v[this.IX(i0,j0)]+t1*v[this.IX(i0,j1)])+
                           s1*(t0*v[this.IX(i1,j0)]+t1*v[this.IX(i1,j1)]));
        
        
                    var speed2=vx*vx+vy*vy; 
                    if(speed2>0.0000001) tnext=.5/Math.sqrt(speed2);
                    else tnext=tleft;
        
                    t=tnext > tleft ? tleft : tnext;
                    tleft-=t;
        
        
                    x+=t*vx;
                    y+=t*vy;
                }
        
        
                if (x<0.5) x=0.5; if (x>N+0.5) x=N+0.5; 
                if (y<0.5) y=0.5; if (y>N+0.5) y=N+0.5; 
                var c = this.checkBounds(x,y);
                i0=c[0]; i1=i0+1;
                j0=c[1]; j1=j0+1;

                s1 = x-i0; s0 = 1-s1; t1 = y-j0; t0 = 1-t1;
                for(var ar = 0 ; ar < d.length; ar ++){
                    d[ar][this.IX(i,j)] =   s0*(t0*d0[ar][this.IX(i0,j0)]+t1*d0[ar][this.IX(i0,j1)])+
                                            s1*(t0*d0[ar][this.IX(i1,j0)]+t1*d0[ar][this.IX(i1,j1)]);            
                }
            }
        }
        //for(var den = 0; den < d.length; den++){
           // this.set_bnd(b, d[den], N);
        //}
    }
    set_bnd(b, x, N){

        if(this.wrap){ 
            
            
            for(var i = 1; i < N - 1; i++) {
                //y axis, top
                var top = this.IX(i, 0  );                
                //y axis bottom
                var bottom = this.IX(i, N-1);
                //x axis left
                var left = this.IX(0 , i);
                //x axis right
                var right = this.IX(N-1, i);            

                x[top]       = x[bottom - N];
                x[bottom]    = x[top + N];
                x[left]      = x[right - 1];
                x[right]     = x[left + 1];
            }
        
        } else {
            for(var i = 1; i < N - 1; i++) {
                ix = this.IX(i, 0  );
                x[ ix ] = b == 2 ? -x[ this.IX(i,j+1)] : x[ this.IX(i,j+1) ];

                ix = this.IX(i, N-1);
                x[ ix ] = b == 2 ? -x[ this.IX(i,j-1) ] : x[ this.IX(i,j-1) ];

                ix = this.IX(0 , i);
                x[ ix ] = b == 1 ? -x[ this.IX(i+1,j) ] : x[ this.IX(i+1,j) ];

                ix = this.IX(N-1, i);            
                x[ ix ] = b == 1 ? -x[ this.IX(i-1,j)] : x[ this.IX(i-1,j)];

            }
        
            
            ix = this.IX(0, 0);
            x[ ix ] = 0.5 * ( x[ this.IX(i+1,j)]      + x[ this.IX(i,j+1) ]       + x[ ix ] );

            ix = this.IX(0, N-1);
            x[ ix ] = 0.5 * ( x[ this.IX(i+1,j) - N ] + x[ this.IX(i,j-1) ]       + x[ ix ]);

            ix = this.IX(N-1, 0);
            x[ ix ] = 0.5 * ( x[ this.IX(i-1,j) ]     + x[ ix -1 + N ]    + x[ ix ]);

            ix = this.IX(N-1, N-1);
            x[ ix ] = 0.5 * ( x[ this.IX(i-1,j) ]     + x[ this.IX(i,j-1) ]       + x[ ix ]);
        }
    }

    edgeVelocities(x,y, N){

        for(var i = 0; i < N - 1; i++) {
            for(var j = 0; j < N - 1; j++) {
                var ix = this.IX(i,j);
                if(!this.solids[ix]){
                    //x
                    //get the direction of force
                    if(x[ix] < 0 ){ //going left, it came from the right, we need to reflect back to the right so cell > 0
                        //first check if the cell to the left is a solid or not, we dont want to reflect velocity into another solid
                        if(this.solids[ix-1]){
                            x[ix] = -x[ix];
                        }
                    }

                    if(x[ix] > 0 ){ //going left, it came from the right, we need to reflect back to the right so cell < 0
                        //first check if the cell to the left is a solid or not, we dont want to reflect velocity into another solid
                        if(this.solids[ix+1]){
                            x[ix] = -x[ix];
                        } 
                    }
                    //y
                    if(y[ix] < 0 ){ //going up, it came from the bottom, we need to reflect back to the top so cell > 0
                        
                        if(this.solids[ix-N]){
                            y[ix] = -y[ix];
                        }
                    }

                    if(y[ix] > 0 ){ //going down, it came from the top, we need to reflect back to down so cell < 0
                        if(this.solids[ix+N]){
                            y[ix] = -y[ix];
                        } 
                    }                   
                }
            }
        }

    }


    addItem(x,y,type,amount){
        var ix = this.IXWrap(x,y);
        if(type == "red"){
            this.densityR[ ix ] += amount;
            //this.densityR[ this.IX(i+1,j) ] += amount;
            //this.densityR[ this.IX(i-1,j) ] += amount;
            //this.densityR[ ix + this.size ] += amount;
            //this.densityR[ ix - this.size ] += amount;
        }
        if(type == "green"){
            this.densityG[ ix ] += amount;
            //this.densityG[ this.IX(i+1,j) ] += amount;
            //this.densityG[ this.IX(i-1,j) ] += amount;
            //this.densityG[ ix + this.size ] += amount;
            //this.densityG[ ix - this.size ] += amount;
        }
        if(type == "blue"){
            this.densityB[ ix ] += amount;
            //this.densityB[ this.IX(i+1,j) ] += amount;
            //this.densityB[ this.IX(i-1,j) ] += amount;
            //this.densityB[ ix + this.size ] += amount;
            //this.densityB[ ix - this.size ] += amount;
        }
        if(type == "addwall"){
            this.solids[ ix ] = true;

        }
        if(type == "removewall"){
            this.solids[ ix ] = false;

        }
    }
    addDensity(x, y, amount){
        this.density[ this.IXWrap(x,y) ] += amount;
        
    }
    addVelocity(x, y, amountX, amountY){
        this.Vx[ this.IXWrap(x,y) ] += amountX;
        this.Vy[ this.IXWrap(x,y) ] += amountY;
    }

    getDensity(x,y){       
        return this.densityR[ this.IXWrap(x,y) ]+','+
                    this.densityG[ this.IXWrap(x,y) ]+','+
                    this.densityB[ this.IXWrap(x,y) ];       
    }

    getVelocityX(x,y){
        
        return this.Vx[ this.IXWrap(x,y) ];
        
    }
    getVelocityY(x,y){
        
        return this.Vy[ this.IXWrap(x,y) ];
        
    }

    IX(x, y){
        return this.ixList[x][y].ix;
    }
    IXObj(x, y){        
        return this.ixList[x][y];
    }
    IXWrap(x, y){
        var c = this.checkBounds(x,y);
        return c[0] + c[1] * this.size;;

    }
    checkBounds(x,y){
        x = Math.floor(x);
        y = Math.floor(y);
        if(x > this.size-1){
            x = x - this.size; 
        } 
        if(x < 0){
            x = x + this.size;
        }

        if(y > this.size-1){
            y = y - this.size; 
        } 
        if(y < 0){
            y = y + this.size;
        }
        return [x,y];
    }
    render(){
        
        this.graphics.clear();
        

        for(var i = 0; i < this.size; i++){
            for(var j = 0; j < this.size; j++){
                var x = i * this.blockSize;
                var y = j * this.blockSize;

                var ix = this.IX(i,j);
                if( this.solids[ ix ] ){
                    var color = Phaser.Display.Color.GetColor(255,255,255);
                    this.graphics.fillStyle(color);
                    this.graphics.fillRect(x, y, this.blockSize, this.blockSize);
                } else {
                    var r = this.densityR[ ix ];
                    if(r > 255){
                        r = 255;
                    }
                    var g = this.densityG[ ix ];
                    if(g > 255){
                        g = 255;
                    }
                    var b = this.densityB[ ix ];
                    if(b > 255){
                        b = 255;
                    }
                    if(r < 10 && b < 10 && g < 10){ continue; }
                    var color = Phaser.Display.Color.GetColor(r,g,b);
                    this.graphics.fillStyle(color);
                    this.graphics.fillRect(x, y, this.blockSize, this.blockSize);
                }
                
            }    
        }
        
        for(var i = 0; i < this.size; i++){
            for(var j = 0; j < this.size; j++){
                
                var x = i ;
                var y = j ;

                var vx = this.Vx[ this.IX(i,j) ] ;
                var vy = this.Vy[ this.IX(i,j) ] ;

                var r = Math.abs(vx * 100);
                var b = Math.abs(vy * 100);
                if(r > 255){
                    r = 255;
                }
                
                if(b > 255){
                    b = 255;
                }
                var color = Phaser.Display.Color.GetColor(r,128,b);
                this.graphics.lineStyle(2,color,0.5);
                //if(!(Math.abs(vx) < 0.1 && Math.abs(vy) <= 0.5)){
                    this.graphics.lineBetween(x* this.blockSize, y* this.blockSize, (x+vx)* this.blockSize,(y+vy)* this.blockSize);
                //} 
                
            }    
        }
        
    }

}

