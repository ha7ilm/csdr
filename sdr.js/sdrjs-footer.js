// ========================================================== 
// ========= / THE CODE COMPILED BY EMCC ENDS HERE ==========
// ========================================================== 

asm$ =
{
	malloc: function(type, size) 
	{
		real_size=size*type.BYTES_PER_ELEMENT;
		pointer = Module._malloc(real_size);
		heap = new Uint8Array(Module.HEAPU8.buffer, pointer, real_size);
		return {
			asm$: true,
			ptr: heap.byteOffset,
			free: function() { Module._free(this.ptr); },
			arr: new type(heap.buffer, heap.byteOffset, size),
			size: size
		};
	},
	cpy: function(dst, dst_offset, src, src_offset, size) 
	{
		if(typeof dst.asm$!='undefined') dst=dst.arr;
		if(typeof src.asm$!='undefined') src=src.arr;
		for(var i=0;i<size;i++) 
			dst[dst_offset+i]=src[src_offset+i];
	}
};

//                                                  void firdes_lowpass_f(float *output, int length, float cutoff_rate, window_t window)
firdes_lowpass_f = Module.cwrap('firdes_lowpass_f', null,                 ['number',     'number',   'number',          'number']);

//                                                            rational_resampler_ff_t rational_resampler_ff(float *input, float *output, int input_size, int interpolation, int decimation, float *taps, int taps_length, int last_taps_delay)
rational_resampler_ff = Module.cwrap('rational_resampler_ff', 'struct',                                    ['number',     'number',      'number',       'number',          'number',       'number',    'number',        'number']);


rational_resampler_ff=function(pinput,poutput,input_length,interpolation,decimation,ptaps,taps_length,last_taps_delay )
{	stackbase=STACKTOP;
	STACKTOP+=4*3;
	_rational_resampler_ff(stackbase, pinput, poutput, input_length, interpolation, decimation, ptaps, taps_length,last_taps_delay);
	returnstruct={ input_processed: getValue(stackbase,'i32'), output_size: getValue(stackbase+4,'i32'), last_taps_delay: getValue(stackbase+8,'i32') };
	STACKTOP=stackbase;
	return returnstruct;
}

sdrjs={};

sdrjs.WINDOW_BOXCAR=0;
sdrjs.WINDOW_BLACKMAN=1;
sdrjs.WINDOW_HAMMING=2;

//this will be impportant whil converting arrays
//http://stackoverflow.com/questions/25839216/convert-float32array-to-int16array

/*sdrjs.prototype.FirdesLowpassF=function(taps_length,transition_bw,window)
{
	this.calculate=function(){}
	this.get_output=function(){}
	this.get_output_heap=function(){}
};*/


sdrjs.ConvertI16_F=function(i16data)
{
	var f32data=new Float32Array(i16data.length);
	for(var i=0;i<i16data.length;i++) f32data[i]=i16data[i]/32768;
	return f32data;
}

ima_adpcm_codec=function(encode,pinput,poutput,input_length,state)
{
	myfunc=(encode)?_encode_ima_adpcm_i16_u8:_decode_ima_adpcm_u8_i16;
	stackbase=STACKTOP;
	STACKTOP+=4*2; //sizeof(int)*2
	myfunc(stackbase, pinput, poutput, input_length, state.ptr);
	state.arr[0]=getValue(stackbase+0,'i32');
	state.arr[1]=getValue(stackbase+4,'i32');
	STACKTOP=stackbase;
};

sdrjs.ImaAdpcm=function()
{
	this.BUFSIZE=1024*64;
	this.ima_adpcm_state=asm$.malloc(Int32Array,2);
	this.i16_buffer=asm$.malloc(Int16Array,this.BUFSIZE*2);
	this.u8_buffer=asm$.malloc(Uint8Array,this.BUFSIZE);
	this.ima_adpcm_state.arr[0]=0;
	this.ima_adpcm_state.arr[1]=0;

	this.encode=function(data)
	{
		//not_tested_yet
		asm$.cpy(this.i16_buffer.arr,0,data,0,data.length);
		ima_adpcm_codec(true,this.i16_buffer.ptr,this.u8_buffer.ptr,data.length,this.ima_adpcm_state);		
		out=new Uint8Array(data.length/2);
		asm$.cpy(out,0,this.u8_buffer,0,data.length/2);
		return out;	
	};

	this.decode=function(data)
	{
		asm$.cpy(this.u8_buffer.arr,0,data,0,data.length);
		ima_adpcm_codec(false,this.u8_buffer.ptr,this.i16_buffer.ptr,data.length,this.ima_adpcm_state);
		out=new Int16Array(data.length*2);
		asm$.cpy(out,0,this.i16_buffer.arr,0,data.length*2);
		return out;	
	};
	this.reset=function() { this.ima_adpcm_state.arr[0]=this.ima_adpcm_state.arr[1]=0|0; }
};

sdrjs.REBUFFER_FIXED=0; //rebuffer should return arrays of fixed size
sdrjs.REBUFFER_MAX=1;	//rebuffer should return arrays with a maximal size of the parameter size

sdrjs.Rebuffer=function(size,mode)
{
	this.mode=mode;
	this.size=size;
	this.total_size=0;
	this.arrays=[];
	this.last_arr=[];
	this.last_arr_offset=0;
	this.push=function(data)
	{
		this.total_size+=data.length;
		this.arrays.push(data);
	};
	this.remaining=function() 
	{ 
		var fixed_bufs_num=Math.floor(this.total_size/this.size);
		if(!this.mode) return fixed_bufs_num;
		else return fixed_bufs_num+(!!(this.total_size-fixed_bufs_num*this.size)); //if REBUFFER_MAX, add one if we could return one more buffer (smaller than the fixed size)
	};
	this.take=function() { var a=this._take(); /*console.log(a);*/ return a; };
	this._take=function() 
	{
		var remain=this.size;
		var offset=0;
		var obuf=new Float32Array(size);
		//console.log("==== get new obuf ====", size);
		while(remain)
		{
			if(this.last_arr_offset==this.last_arr.length)
			{
				if(this.arrays.length==0)
				{ 
					//console.log("this should not happen");
					if(this.mode) //REBUFFER_MAX
					{
						this.total_size=0;
						return obuf.subarray(0,offset);  
					}
					else return new Float32Array(0); //REBUFFER_FIXED
				}
				//console.log("pick new last_arr");
				this.last_arr=this.arrays.shift();
				this.last_arr_offset=0;
			}
			var rwithin=this.last_arr.length-this.last_arr_offset;
			//console.log("b :: ","remain", remain, "rwithin",rwithin,"last_arr.length",this.last_arr.length,"larroffset",this.last_arr_offset,"offset",offset);
			if(remain<rwithin)
			{
				//console.log("remain < rwithin"); //seems problematic @Andris
				for(var i=0;i<remain;i++) obuf[offset++]=this.last_arr[this.last_arr_offset++];
				remain=0;
			}
			else
			{
				//console.log("remain > rwithin");
				for(var i=0;i<rwithin;i++) obuf[offset++]=this.last_arr[this.last_arr_offset++];
				remain-=rwithin;
			}
			//console.log("e :: ","remain", remain, "rwithin",rwithin,"last_arr.length",this.last_arr.length,"larroffset",this.last_arr_offset,"offset",offset);
		}
			
		this.total_size-=obuf.length;
		//console.log("return _take");
		return obuf;
	};
};

sdrjs.RationalResamplerFF=function(interpolation,decimation,transition_bw,window)
{
	this.interpolation=interpolation;
	this.decimation=decimation;
	this.transition_bw = (typeof transition_bw=='undefined')?0.05:transition_bw;
	this.window = (typeof window=='undefined')?1:window;
	this.buffer_size=1024*512;
	this.output_buffer_size=Math.floor((this.buffer_size*interpolation)/decimation);
	this.input_buffer = asm$.malloc(Float32Array,this.buffer_size);
	this.output_buffer = asm$.malloc(Float32Array,this.output_buffer_size);
	//Calculate filter
	this.taps_length = Math.floor(4/this.transition_bw);
	this.taps = asm$.malloc(Float32Array,this.taps_length);
	var cutoff_for_interpolation=1.0/interpolation;
	var cutoff_for_decimation=1.0/decimation;
	var cutoff = (cutoff_for_interpolation<cutoff_for_decimation)?cutoff_for_interpolation:cutoff_for_decimation; //get the lower
	firdes_lowpass_f(this.taps.ptr, this.taps_length, cutoff/2, window);
	
	this.remain = 0;
	this.remain_offset=0;
	this.last_taps_delay=0;

	this.process=function(input)
	{
		
		if(input.length+this.remain > this.buffer_size)
		{
			return new Float32Array(0); console.log("sdrjs.RationalResamplerFF: critical audio buffering error"); //This should not happen...
		/*	console.log("RationalResamplerFF: splitting..."); //TODO: this branch has not been checked
			output_buffers=Array();
			new_buffer_size=this.buffer_size/2;
			i=0;
			//process the input in chunks of new_buffer_size, and add the output product Float32Array-s to output_buffers.
			while((i++)*new_buffer_size<=input.length)
			{
				output_buffers.push(this._process_noheapcheck(input.subarray(i*new_buffer_size,(i+1)*new_buffer_size)));
			}
			//add up the sizes of the output_buffer-s.
			total_output_length=0;
			output_buffers.forEach(function(a){total_output_length+=a.length;});
			//create one big buffer from concatenating the output_buffer-s
			output=new Float32Array(total_output_length);
			output_pos=0;
			output_buffers.forEach(function(a){
				asm$.cpy(output,output_pos,a,0,a.length);
				output_pos+=a.length;
			});
			return output;*/
		}
		else return this._process_noheapcheck(input);
	};
	this._process_noheapcheck=function(input) //if we are sure we have enough space in the buffers 
	{
		asm$.cpy(this.input_buffer.arr,0,this.input_buffer.arr,this.remain_offset,this.remain);
		asm$.cpy(this.input_buffer.arr, this.remain, input, 0, input.length);
		var total_input_size=input.length+this.remain;
		d=rational_resampler_ff(this.input_buffer.ptr, this.output_buffer.ptr, total_input_size, this.interpolation, this.decimation, this.taps.ptr, this.taps_length, this.last_taps_delay);
		this.last_taps_delay=d.last_taps_delay;
		this.remain=total_input_size-d.input_processed;
		this.remain_offset=d.input_processed;
		var output_copy_arr=new Float32Array(d.output_size);
		asm$.cpy(output_copy_arr,0,this.output_buffer.arr,0,d.output_size);
		return output_copy_arr;
	};
};


_sdrjs_logb=function(what) { document.body.innerHTML+=what+"<br />"; }


function test_firdes_lowpass_f_original()
{
	//Original method explained over here: 
	//http://kapadia.github.io/emscripten/2013/09/13/emscripten-pointers-and-pointers.html
	_sdrjs_logb("test_firdes_lowpass_f_original():");
	_sdrjs_logb("Now designing FIR filter with firdes_lowpass_f in sdr.js...");
	_sdrjs_logb("output should be the same as: <strong>csdr firdes_lowpass_f 0.1 101 HAMMING</strong>");
	
	var outputSize = 101*4;
	var outputPtr = Module._malloc(outputSize);
	var outputHeap = new Uint8Array(Module.HEAPU8.buffer, outputPtr, outputSize);
	firdes_lowpass_f(outputHeap.byteOffset,101,0.1,2);
	var output = new Float32Array(outputHeap.buffer, outputHeap.byteOffset, 101);
	outputStr=String();
	for(i=0;i<output.length;i++) outputStr+=output[i].toFixed(6)+", ";
	Module._free(outputHeap.byteOffset);
	_sdrjs_logb(outputStr);
}


function test_firdes_lowpass_f_new()
{
	//This is much simpler, using asm$
	_sdrjs_logb("test_firdes_lowpass_f_new():");
	_sdrjs_logb("Now designing FIR filter with firdes_lowpass_f in sdr.js...");
	_sdrjs_logb("output should be the same as: <strong>csdr firdes_lowpass_f 0.1 101 HAMMING</strong>");
	
	output=asm$.malloc(Float32Array,101);
	firdes_lowpass_f(output.ptr,101,0.1,2);
	outputStr=String();
	for(i=0;i<output.arr.length;i++) outputStr+=(output.arr[i]).toFixed(6)+", ";
	output.free();
	_sdrjs_logb(outputStr);
}

function test_struct_return_value()
{
	v=STACKTOP;
	STACKTOP+=4*3;
	_shift_addition_init(v,0.2);
	console.log( 
		"sinval=", getValue(v,'float'), 
		"cosval=", getValue(v+4,'float'), 
		"rate=", getValue(v+8,'float') 
	);
	STACKTOP=v;
}
