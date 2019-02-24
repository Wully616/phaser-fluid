//https://mikeash.com/pyblog/fluid-simulation-for-dummies.html
//http://www.cs.northwestern.edu/~sco590/cs140/assignment3.html
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
        this.createDebugGrid(this.scene,this.gridWidth,this.gridHeight,this.blockSize);
        
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


        this.solids = new Array(this.gridsize).fill(1);

        this.graphics = this.scene.add.graphics();

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
            
        this.diffuseVec(1, Vx0, Vx, 2, Vy0, Vy, visc, dt, this.iter, N );
        
        this.project(Vx0, Vy0, Vx, Vy, this.iter, N);
        
        this.advectVec(1, Vx, Vx0, 2, Vy, Vy0, Vx0, Vy0,  dt, N); 
        
        this.project(Vx, Vy, Vx0, Vy0, this.iter, N);
        
        //process R G B channels of density
        this.diffuseRGB(0, [sR,sG,sB], [densityR,densityG,densityB], diff, dt, this.iter, N);
        this.advectRGB(0, [densityR,densityG,densityB], [sR,sG,sB], Vx, Vy, dt, N);

  
    }

    //type, oldvel, newvel, visc, timestep, iteration, size
    diffuse(b, x, x0, diff, dt, iter, N){
        var a = dt * diff * (N - 2) * (N - 2);
        this.lin_solve(b, x, x0, a, 1 + 4 * a, iter, N);
    }
    diffuseVec(bndsA, curA, prevA, bndsB, curB, prevB, diff, dt, iter, N){
        var a = dt * diff * (N - 2) * (N - 2);
        this.lin_solveVec(bndsA,curA, prevA, bndsB, curB, prevB, a, 1 + 6 * a, iter, N);
    }
    diffuseRGB(b, x, x0, diff, dt, iter, N){
        var a = dt * diff * (N - 2) * (N - 2);
        this.lin_solveRGB(b, x, x0, a, 1 + 4 * a, iter, N);
    }

    lin_solve(b, x, x0,  a,  c,  iter,  N){
        var cRecip = 1.0 / c;
        
        for (var k = 0; k < iter; k++) {     
            for (var j = 1; j < N - 1; j++) {
                for (var i = 1; i < N - 1; i++) {
                    var ix = this.IX(i, j);
                    x[ ix ] = this.calcNeighborVal(ix ,x,x0,a, cRecip, N);                    
                }
            }            
            this.set_bnd(b, x, N);
        }
    }
    lin_solveRGB(b, x, x0,  a,  c,  iter,  N){
        var cRecip = 1.0 / c;
        for (var k = 0; k < iter; k++) {     
            for (var j = 1; j < N - 1; j++) {
                for (var i = 1; i < N - 1; i++) {
                    var ix = this.IX(i, j);
                    x[0][ ix ] = this.calcNeighborVal(ix ,x[0],x0[0],a, cRecip, N); //r
                    x[1][ ix ] = this.calcNeighborVal(ix ,x[1],x0[1],a, cRecip, N); //g
                    x[2][ ix ] = this.calcNeighborVal(ix ,x[2],x0[2],a, cRecip, N); //b

                }
            }      
            
            this.set_bnd(b, x[0], N);
            this.set_bnd(b, x[1], N);
            this.set_bnd(b, x[2], N);
            
        }
    }
    lin_solveVec(b, x, x0, b2, y, y0,  a,  c,  iter,  N){
        var cRecip = 1.0 / c;
        for (var k = 0; k < iter; k++) {     
            for (var j = 1; j < N - 1; j++) {
                for (var i = 1; i < N - 1; i++) {
                    var ix = this.IX(i, j);
                    x[ ix ] = this.calcNeighborVal(ix ,x,x0,a, cRecip, N);     
                    y[ ix ] = this.calcNeighborVal(ix ,y,y0,a, cRecip, N);                                           
                }
            }            
            this.set_bndVec(b, x, b2, y, N)
        }
    }

    calcNeighborVal(ix ,x,x0,a, cRecip, N){
        return (
                x0[ix] + a *
                (   x[ ix + 1 ]  + //+ 1 on x axis
                    x[ ix - 1 ]  + //-1 on x axis
                    x[ ix + N]  + /// +1 on y axis
                    x[ ix - N] // -1 on y axis
                )
            ) * cRecip;
    }

    project(velocX, velocY, p, div, iter, N)
    {
        var nRecip = 1 / N;   
        for (var j = 1; j < N - 1; j++) {
            for (var i = 1; i < N - 1; i++) {
                var ix = this.IX(i, j);
                div[ ix ] = -0.5*(
                            (velocX[ ix + 1] )
                            -(velocX[ ix -1 ] )
                            +(velocY[ ix + N ])
                            -(velocY[ ix - N ] )
                ) * nRecip;
                p[ ix ] = 0;
            }
        }
    
        this.set_bndVec(0, div, 0, p, N)
        this.lin_solve(0, p, div, 1, 6, iter, N);
        
        
        for (var j = 1; j < N - 1; j++) {
            for (var i = 1; i < N - 1; i++) {
                var ix = this.IX(i, j);
                velocX[ ix ] -= 0.5 * (  p[ ix + 1 ]- p[ ix - 1 ]* this.solids[ix - 1]) * N;
                velocY[ ix ] -= 0.5 * (  p[ ix + N ]- p[ ix - N ]* this.solids[ix - N]) * N;
            }
        }
        
        this.set_bndVec(1, velocX, 2, velocX, N)
    }

    //refactor
    advectVec(b, d, d0, b2, e, e0, velocX, velocY, dt, N)
    {
        var i0, i1, j0, j1;
        
        var dtx = dt * (N - 2);
        var dty = dt * (N - 2);
        
        var s0, s1, t0, t1;
        var tmp1, tmp2, x, y;
        
        var Nfloat = N;
        var ifloat, jfloat;
        var i, j;
        

        for(j = 1, jfloat = 1; j < N - 1; j++, jfloat++) { 
            for(i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
                tmp1 = dtx * velocX[this.IX(i, j)];
                tmp2 = dty * velocY[this.IX(i, j)];

                x    = ifloat - tmp1; 
                y    = jfloat - tmp2;
                
                if(x < 0.5) x = 0.5; 
                if(x > Nfloat + 0.5) x = Nfloat + 0.5; 
                i0 = Math.floor(x); 
                i1 = i0 + 1.0;
                if(y < 0.5) y = 0.5; 
                if(y > Nfloat + 0.5) y = Nfloat + 0.5; 
                j0 = Math.floor(y);
                j1 = j0 + 1.0; 
                
                s1 = x - i0; 
                s0 = 1.0 - s1; 
                t1 = y - j0; 
                t0 = 1.0 - t1;
                
                var i0i = Math.floor(i0);
                var i1i = Math.floor(i1);
                var j0i = Math.floor(j0);
                var j1i = Math.floor(j1);
                
                var ix1 = this.IX(i0i, j0i);
                var ix2 = this.IX(i0i, j1i);
                var ix3 = this.IX(i1i, j0i);
                var ix4 = this.IX(i1i, j1i);

                d[this.IX(i, j)] = 
                
                    s0 * ( t0 * d0[ix1]  +  t1 * d0[ix2]  ) +
                    s1 * ( t0 * d0[ix3]  +  t1 * d0[ix4] );
                
                e[this.IX(i, j)] = 
                
                    s0 * ( t0 * e0[ix1]   +  t1 * e0[ix2]  ) +
                    s1 * ( t0 * e0[ix3]  +  t1 * e0[ix4] );                    
            }
        }
    
        this.set_bndVec(b, d, b2, e, N)
    }

    //refactor
    advect(b, d, d0,  velocX, velocY, dt, N)
    {
        var i0, i1, j0, j1;
        
        var dtx = dt * (N - 2);
        var dty = dt * (N - 2);
        
        var s0, s1, t0, t1;
        var tmp1, tmp2, x, y;
        
        var Nfloat = N;
        var ifloat, jfloat;
        var i, j;
        

        for(j = 1, jfloat = 1; j < N - 1; j++, jfloat++) { 
            for(i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
                tmp1 = dtx * velocX[this.IX(i, j)];
                tmp2 = dty * velocY[this.IX(i, j)];

                x    = ifloat - tmp1; 
                y    = jfloat - tmp2;
                
                if(x < 0.5) x = 0.5; 
                if(x > Nfloat + 0.5) x = Nfloat + 0.5; 
                i0 = Math.floor(x); 
                i1 = i0 + 1.0;
                if(y < 0.5) y = 0.5; 
                if(y > Nfloat + 0.5) y = Nfloat + 0.5; 
                j0 = Math.floor(y);
                j1 = j0 + 1.0; 
                
                s1 = x - i0; 
                s0 = 1.0 - s1; 
                t1 = y - j0; 
                t0 = 1.0 - t1;
                
                var i0i = Math.floor(i0);
                var i1i = Math.floor(i1);
                var j0i = Math.floor(j0);
                var j1i = Math.floor(j1);
                var ix1 = this.IX(i0i, j0i);
                var ix2 = this.IX(i0i, j1i);
                var ix3 = this.IX(i1i, j0i);
                var ix4 = this.IX(i1i, j1i);
                d[this.IX(i, j)] = 
                
                    s0 * ( t0 * d0[ix1]  +  t1 * d0[ix2]  ) +
                    s1 * ( t0 * d0[ix3]  +  t1 * d0[ix4] );
            }
        }
    
        this.set_bnd(b, d, N);
    }
    advectRGB(b, d, d0,  velocX, velocY, dt, N)
    {
        var i0, i1, j0, j1;
        
        var dtx = dt * (N - 2);
        var dty = dt * (N - 2);
        
        var s0, s1, t0, t1;
        var tmp1, tmp2, x, y;
        
        var Nfloat = N;
        var ifloat, jfloat;
        var i, j;
        

        for(j = 1, jfloat = 1; j < N - 1; j++, jfloat++) { 
            for(i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
                tmp1 = dtx * velocX[this.IX(i, j)];
                tmp2 = dty * velocY[this.IX(i, j)];

                x    = ifloat - tmp1; 
                y    = jfloat - tmp2;
                
                if(x < 0.5) x = 0.5; 
                if(x > Nfloat + 0.5) x = Nfloat + 0.5; 
                i0 = Math.floor(x); 
                i1 = i0 + 1.0;
                if(y < 0.5) y = 0.5; 
                if(y > Nfloat + 0.5) y = Nfloat + 0.5; 
                j0 = Math.floor(y);
                j1 = j0 + 1.0; 
                
                s1 = x - i0; 
                s0 = 1.0 - s1; 
                t1 = y - j0; 
                t0 = 1.0 - t1;
                
                var i0i = Math.floor(i0);
                var i1i = Math.floor(i1);
                var j0i = Math.floor(j0);
                var j1i = Math.floor(j1);
                var ix1 = this.IX(i0i, j0i);
                var ix2 = this.IX(i0i, j1i);
                var ix3 = this.IX(i1i, j0i);
                var ix4 = this.IX(i1i, j1i);
                for(var den = 0; den < d.length; den++){
                    d[den][this.IX(i, j)] = 
                
                    s0 * ( t0 * d0[den][ix1] +  t1 * d0[den][ix2] ) +
                    s1 * ( t0 * d0[den][ix3] +  t1 * d0[den][ix4] );

                }
            }
        }
        for(var den = 0; den < d.length; den++){
            this.set_bnd(b, d[den], N);
        }
    }
    set_bndVec(b, x, b2, y, N){

        var ix;
        for(var i = 1; i < N - 1; i++) {
            //top edge
            ix = this.IX(i, 0  );
            x[ ix ] = b == 2 ? -x[ ix + N] : x[ ix + N ];
            y[ ix ] = b2 == 2 ? -y[ ix + N] : y[ ix + N ];
            //bottom edge
            ix = this.IX(i, N-1);
            x[ ix ] = b == 2 ? -x[ ix - N ] : x[ ix - N ];
            y[ ix ] = b2 == 2 ? -y[ ix - N ] : y[ ix - N ];
            //left edge
            ix = this.IX(0 , i);
            x[ ix ] = b == 1 ? -x[ ix + 1 ] : x[ ix + 1 ];
            y[ ix ] = b2 == 1 ? -y[ ix + 1 ] : y[ ix + 1 ];
            //right edge
            ix = this.IX(N-1, i);            
            x[ ix ] = b == 1 ? -x[ ix - 1] : x[ ix - 1];
            y[ ix ] = b2 == 1 ? -y[ ix - 1] : y[ ix - 1];
 
        }

        
        ix = this.IX(0, 0);
        x[ ix ] = 0.5 * ( x[ ix + 1]      + x[ ix + N ]       + x[ ix ] );
        y[ ix ] = 0.5 * ( y[ ix + 1]      + y[ ix + N ]       + y[ ix ] );

        ix = this.IX(0, N-1);
        x[ ix ] = 0.5 * ( x[ ix + 1 - N ] + x[ ix - N ]       + x[ ix ]);
        y[ ix ] = 0.5 * ( y[ ix + 1 - N ] + y[ ix - N ]       + y[ ix ]);

        ix = this.IX(N-1, 0);
        x[ ix ] = 0.5 * ( x[ ix - 1 ]     + x[ ix -1 + N ]    + x[ ix ]);
        y[ ix ] = 0.5 * ( y[ ix - 1 ]     + y[ ix -1 + N ]    + y[ ix ]);

        ix = this.IX(N-1, N-1);
        x[ ix ] = 0.5 * ( x[ ix - 1 ]     + x[ ix - N ]       + x[ ix ]);
        y[ ix ] = 0.5 * ( y[ ix - 1 ]     + y[ ix - N ]       + y[ ix ]);
    }

    set_bnd(b, x, N){

        var ix;
        for(var i = 1; i < N - 1; i++) {
            ix = this.IX(i, 0  );
            x[ ix ] = b == 2 ? -x[ ix + N] : x[ ix + N ];

            ix = this.IX(i, N-1);
            x[ ix ] = b == 2 ? -x[ ix - N ] : x[ ix - N ];

            ix = this.IX(0 , i);
            x[ ix ] = b == 1 ? -x[ ix + 1 ] : x[ ix + 1 ];

            ix = this.IX(N-1, i);            
            x[ ix ] = b == 1 ? -x[ ix - 1] : x[ ix - 1];

        }
    
        
        ix = this.IX(0, 0);
        x[ ix ] = 0.5 * ( x[ ix + 1]      + x[ ix + N ]       + x[ ix ] );

        ix = this.IX(0, N-1);
        x[ ix ] = 0.5 * ( x[ ix + 1 - N ] + x[ ix - N ]       + x[ ix ]);

        ix = this.IX(N-1, 0);
        x[ ix ] = 0.5 * ( x[ ix - 1 ]     + x[ ix -1 + N ]    + x[ ix ]);

        ix = this.IX(N-1, N-1);
        x[ ix ] = 0.5 * ( x[ ix - 1 ]     + x[ ix - N ]       + x[ ix ]);
    }

    addItem(x,y,type,amount){
        var ix = this.IX(x,y);
        if(type == "red"){
            this.densityR[ ix ] += amount;
            this.densityR[ ix + 1 ] += amount;
            this.densityR[ ix - 1 ] += amount;
            this.densityR[ ix + this.size ] += amount;
            this.densityR[ ix - this.size ] += amount;
        }
        if(type == "green"){
            this.densityG[ ix ] += amount;
            this.densityG[ ix + 1 ] += amount;
            this.densityG[ ix - 1 ] += amount;
            this.densityG[ ix + this.size ] += amount;
            this.densityG[ ix - this.size ] += amount;
        }
        if(type == "blue"){
            this.densityB[ ix ] += amount;
            this.densityB[ ix + 1 ] += amount;
            this.densityB[ ix - 1 ] += amount;
            this.densityB[ ix + this.size ] += amount;
            this.densityB[ ix - this.size ] += amount;
        }
        if(type == "addwall"){
            this.solids[ ix ] = 0;
            this.solids[ ix + 1 ] = 0;
            this.solids[ ix - 1 ] = 0;
            this.solids[ ix + this.size ] = 0;
            this.solids[ ix - this.size ] = 0;
            this.solids[ ix + 1 + this.size] = 0;
            this.solids[ ix + 1 - this.size] = 0;
            this.solids[ ix - 1 + this.size ] = 0;
            this.solids[ ix - 1 - this.size ] = 0;
        }
        if(type == "removewall"){
            this.solids[ ix ] = 1;
            this.solids[ ix + 1 ] = 1;
            this.solids[ ix - 1 ] = 1;
            this.solids[ ix + this.size ] = 1;
            this.solids[ ix - this.size ] = 1;
            this.solids[ ix + 1 + this.size] = 1;
            this.solids[ ix + 1 - this.size] = 1;
            this.solids[ ix - 1 + this.size ] = 1;
            this.solids[ ix - 1 - this.size ] = 1;
        }
    }
    addDensity(x, y, amount){
        this.density[ this.IX(x,y) ] += amount;
        
    }
    addVelocity(x, y, amountX, amountY){
        this.Vx[ this.IX(x,y) ] += amountX;
        this.Vy[ this.IX(x,y) ] += amountY;
    }

    getDensity(x,y){
        return this.densityR[ this.IX(x,y) ]+','+
                this.densityG[ this.IX(x,y) ]+','+
                this.densityB[ this.IX(x,y) ];
    }


    getVelocityX(x,y){
        return this.Vx[ this.IX(x,y) ];
    }
    getVelocityY(x,y){
        return this.Vy[ this.IX(x,y) ];
    }

    IX(x, y){
        //converts the X, Y coords to a 1D array, size being the length of Y

        return Math.max(0, Math.min( x, this.size-1)) + Math.max(0, Math.min( y, this.size-1)) * this.size;
    }

    
    render(){
        
        this.graphics.clear();
        

        for(var i = 0; i < this.size; i++){
            for(var j = 0; j < this.size; j++){
                var x = i * this.blockSize;
                var y = j * this.blockSize;

                var ix = this.IX(i,j);
                if( this.solids[ ix ] == 0 ){
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

                var vx = this.Vx[ this.IX(i,j) ];
                var vy = this.Vy[ this.IX(i,j) ];

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
                if(!(Math.abs(vx) < 0.1 && Math.abs(vy) <= 0.5)){
                    this.graphics.lineBetween(x* this.blockSize, y* this.blockSize, (x+vx)* this.blockSize,(y+vy)* this.blockSize);
                } 
                
            }    
        }
        
    }

}

